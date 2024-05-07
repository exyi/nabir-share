#!/usr/bin/env python3

import functools
import itertools
from multiprocessing.pool import Pool
from typing import Any, Callable, Generator, Iterator, List, Literal, Optional, Sequence, Tuple, TypeVar, Union
import Bio.PDB, Bio.PDB.Structure, Bio.PDB.Model, Bio.PDB.Residue, Bio.PDB.Atom, Bio.PDB.Chain, Bio.PDB.Entity
import os, sys, io, gzip, math, json, numpy as np
import polars as pl
from dataclasses import dataclass
import dataclasses
from requests import HTTPError
import multiprocessing
import scipy.spatial
from scipy.spatial.transform import Rotation
from pair_csv_parse import scan_pair_csvs
import pdb_utils
import pair_defs as pair_defs
from async_utils import MockPool
pdef = pair_defs

_sentinel = object()

@dataclass
class AltResidue:
    """
    Residue wrapper filtering only atoms with a given altloc code
    """
    res: Bio.PDB.Residue.Residue
    alt: str

    @property
    def resname(self):
        return self.res.resname

    def transform_atom(self, a: Bio.PDB.Atom.DisorderedAtom) -> Any:
        """unwrap a disordered atom"""
        if a.is_disordered():
            if self.alt == '':
                # don't crash on uninteresting atoms
                if a.id.startswith('H') or a.id in ['P', 'OP1', 'OP2', 'OP3', "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "O2'", "C2'", "C1'"]:
                    return a.disordered_get(a.disordered_get_id_list()[0])
                elif np.allclose(a.disordered_get('A').coord, a.disordered_get('B').coord, atol=0.05):
                    return a.disordered_get('A')
            if not a.disordered_has_id(self.alt):
                raise KeyError(f"Atom {self.res.full_id}:{a.id} is disordered, but alt='{self.alt}' must be one of {a.disordered_get_id_list()}")
            return a.disordered_get(self.alt)
        else:
            return a
    def get_atoms(self):
        """List all residue atoms"""
        disorder_match = False
        for a in self.res.get_atoms():
            if self.alt != '' and a.is_disordered():
                if a.disordered_has_id(self.alt):
                    disorder_match = True
                    yield a.disordered_get(self.alt)
                    continue

            yield self.transform_atom(a)

        if self.alt != '' and not disorder_match:
            raise KeyError(f"Alt {self.alt} not found in residue {self.res.full_id}")
    def get_atom(self, name, default: Any=_sentinel) -> Bio.PDB.Atom.Atom:
        """Get residue atom by name"""
        if name in self.res:
            return self.transform_atom(self.res[name])
        else:
            if default is _sentinel:
                raise KeyError(f"Atom {name} not found in residue {self.res.full_id}")
            else:
                return default
            
    def copy(self) -> 'AltResidue':
        return AltResidue(self.res.copy(), self.alt)

def try_get_base_atoms(res: AltResidue) -> Optional[Tuple[Bio.PDB.Atom.Atom, List[Bio.PDB.Atom.Atom]]]:
    """
    Get atoms af the nucleobase, excluding sugar and phosphate atoms
    return tuple of (C1', [base atoms]), or None in case of an error (C1' doesn't exist)
    """
    c1 = res.get_atom("C1'", None)
    if c1 is not None:
        planar_base_atoms = [ a for a in res.get_atoms() if not a.name.endswith("'") and a.element not in [ "H", "D", "P" ] and a.name not in ["P", "OP1", "OP2", "OP3"] ]
        assert len(planar_base_atoms) >= 6
        return c1, planar_base_atoms
    return None

def get_base_atoms(res: AltResidue) -> Tuple[Bio.PDB.Atom.Atom, List[Bio.PDB.Atom.Atom]]:
    """
    Get atoms af the nucleobase, excluding sugar and phosphate atoms
    return tuple of (C1', [base atoms])
    """
    x = try_get_base_atoms(res)
    if x is None:
        raise ResideTransformError(f"No C1' atom found at {res.res.full_id}")
    return x

def bfs(adj_matrix, current):
    assert np.diagonal(adj_matrix).all()
    while True:
        next = adj_matrix @ current
        if np.array_equal(next, current):
            return next
        current = next

class ResideTransformError(Exception):
    pass

def filter_atoms(c1, atoms):
    """
    Returns only atoms connected to the `c1` atom
    """
    arr = np.array([ a.coord for a in atoms ]) - c1.coord
    c1_bonded = np.linalg.norm(arr, axis=1) <= 1.8
    if not np.any(c1_bonded):
        raise ResideTransformError(f"No bonding atom found at {c1.full_id}")
    adj_matrix = np.linalg.norm(arr[:, None, :] - arr[None, :, :], axis=2) <= 1.6
    connected_component = bfs(adj_matrix, c1_bonded)

    filter = connected_component
    return [ a for a, b in zip(atoms, filter) if b ]

@dataclass
class TranslationThenRotation:
    translation: np.ndarray
    rotation: np.ndarray
    def __post_init__(self):
        assert self.rotation.shape == (3, 3)
        assert self.translation.shape == (3,)
        assert np.isclose(np.linalg.det(self.rotation), 1, atol=1e-4)


@dataclass
class TranslationRotationTranslation:
    translation1: np.ndarray
    rotation: np.ndarray
    translation2: np.ndarray
    def __post_init__(self):
        assert self.rotation.shape == (3, 3)
        assert self.translation1.shape == (3,)
        assert self.translation2.shape == (3,)
        assert np.isclose(np.linalg.det(self.rotation), 1, atol=1e-4)

@dataclass
class ResiduePosition:
    # translation of C1 to origin
    origin: np.ndarray
    # rotation around C1 to align C1'-N bond to x-axis and the plane to X/Y axes
    rotation: np.ndarray
    # orthonormal basis of the base plane (fitted through least squares)
    plane_basis: np.ndarray
    # normal vector of the base plane
    plane_normal_vector: np.ndarray
    # projection matrix to the plane
    projection_matrix: np.ndarray
    # C1' atom coordinates
    c1_pos: Optional[np.ndarray] = None
    # N1/N9 atom coordinates
    n_pos: Optional[np.ndarray] = None

    def __post_init__(self):
        assert self.rotation.shape == (3, 3)
        assert self.plane_basis.shape == (3, 2)

    def plane_projection(self, point: np.ndarray) -> np.ndarray:
        """
        Projects the point onto the base plane
        """
        assert point.shape[-1] == 3
        tpoint = point - self.origin
        return np.matmul(tpoint, self.projection_matrix) + self.origin

    def plane_get_deviation(self, point: np.ndarray) -> float:
        """
        Returns the distance of the point from the plane (in Å), negative if it's on opposite direction of the normal vector
        """
        assert point.shape == (3,)
        distance = _linalg_norm(point - self.plane_projection(point))
        sign = np.sign(np.dot(point - self.origin, self.plane_normal_vector))
        return sign * distance
    
    def get_projection_squish_angle(self, point1: np.ndarray, point2: np.ndarray) -> float:
        """
        Returns the ratio of the distance between the projections and the distance between the points.
        """
        assert point1.shape == (3,) == point2.shape
        proj1 = self.plane_projection(point1)
        proj2 = self.plane_projection(point2)
        dist_cmp = _linalg_norm(proj1 - proj2) / _linalg_norm(point1 - point2)
        return min(1, max(0, dist_cmp))
    
    def get_relative_line_rotation(self, point1: np.ndarray, point2: np.ndarray) -> float:
        """
        Returns the angle between the line and the plane.
        Negative if the line points in the opposite direction than the normal vector.
        """
        assert point1.shape == (3,) == point2.shape
        line = point2 - point1
        normal_angle = math.degrees(math.acos(np.dot(_linalg_normalize(line), self.plane_normal_vector)))
        return 90 - normal_angle

def fit_plane_magic_svd(atoms: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """
    Returns centroid + S (= basis1, basis2, normal)
    https://math.stackexchange.com/a/99317
    """
    assert atoms.shape[1] == 3
    left_singular_vectors, _, _ = np.linalg.svd((atoms - np.mean(atoms, axis=0, keepdims=True)).T)
    assert left_singular_vectors.shape == (3, 3)
    return np.mean(atoms, axis=0), left_singular_vectors

def get_residue_posinfo(res: AltResidue) -> ResiduePosition:
    """
    Returns a tuple of (translation, rotation) that moves the C1' atom to the origin and the C1'-N bond to the x-axis and all planar atoms with (near) 0 z-coordinate.
    """
    c1, planar_atoms = get_base_atoms(res)
    nX = res.get_atom("N1", None) or res.get_atom("N9", None)
    planar_atoms = filter_atoms(c1, planar_atoms)
    
    atom_names = np.array([ a.name for a in planar_atoms ])
    atom_elements = np.array([ a.element for a in planar_atoms ])
    atoms = np.array([ a.coord for a in planar_atoms ])

    centroid, rot_matrix1 = fit_plane_magic_svd(atoms)
    plane_basis = rot_matrix1[:, :2]
    normal_v = rot_matrix1[:, 2]
    atoms_c = atoms - centroid
    sum_error = np.sum(np.dot(atoms_c, normal_v) ** 2)
    if sum_error > 0.5 * len(planar_atoms): # TODO: reasonable threshold?
        raise ResideTransformError(f"Residue really doesn't seem planar: {res.res.full_id}, plane RMSE = {np.sqrt(sum_error / len(planar_atoms))}, planar atoms = {list(atom_names)}")
    
    # project N onto the plane, this will be the origin
    assert plane_basis.shape == (3, 2)
    projection = plane_basis @ plane_basis.T
    if nX is None:
        origin = centroid
        print(f"WARNING: No N1/N9 atom found in {res.res.full_id}.")
        rot = rot_matrix1
    else:
        origin = ((nX.coord - centroid) @ projection) + centroid


        # rotate to align the fitted_plane to x, y
        # rotation matrix is an orthonormal matrix with det=1. Cross product gives us the remaining vector orthogonal to the plane.
        assert np.sum((atoms_c @ rot_matrix1)[:, 2] ** 2) < 1.3 * sum_error
        
        assert np.allclose(rot_matrix1 @ rot_matrix1.T, np.eye(3), atol=1e-5)
        assert np.allclose(np.linalg.det(rot_matrix1), 1, atol=1e-5)

        atoms_c_rot1 = atoms_c @ rot_matrix1

        assert np.allclose(np.linalg.norm(atoms_c, axis=1), np.linalg.norm(atoms_c_rot1, axis=1), atol=1e-5)

        # orient C1'-N bond to x-axis (using 2D rotation on X,Y)
        c1_c_rot = (c1.coord - origin) @ rot_matrix1
        x_vector = -_linalg_normalize(c1_c_rot)
        rot_matrix2 = np.array([
            [x_vector[0], -x_vector[1], 0],
            [x_vector[1], x_vector[0],  0],
            [ 0,          0,            1]
        ])
        rot = rot_matrix1 @ rot_matrix2

        left_c = res.get_atom("C8", None) if nX.name == "N9" else res.get_atom("C6", None)
        if left_c is not None and ((left_c.coord - origin) @ rot)[1] > 0:
            # flip the Y axis
            rot = rot @ np.array([ [1, 0, 0], [0, -1, 0], [0, 0, 1] ])
            normal_v = -normal_v


    return ResiduePosition(
        origin,
        rot,
        plane_basis,
        normal_v,
        projection,
        c1.coord,
        nX.coord if nX is not None else None
    )

def get_residue_posinfo_C1_N(res: AltResidue) -> TranslationThenRotation:
    """
    Orients residue so that the C1-N bond is aligned with the Y-axis, the N1-C2 or N9-C8 bond is towards X axis and remaining is Z
    """
    c1 = res.get_atom("C1'", None)
    if res.get_atom("N9", None) is not None:
        n = res.get_atom("N9", None)
        c2 = res.get_atom("C8", None)
    else:
        n = res.get_atom("N1", None)
        c2 = res.get_atom("C6", None)
    if c1 is None or n is None or c2 is None:
        raise ResideTransformError(f"Missing atoms in residue {res.res.full_id}")
    translation = -n.coord # N1/N9 is origin
    # y axis is aligned with C1'-N
    y = n.coord - c1.coord
    y /= np.linalg.norm(y)
    # x axis is aligned with N1-C2/N9-C8 (but perpendicular to y)
    x = n.coord - c2.coord
    x -= np.dot(y, x) * y # project X onto Y
    x /= np.linalg.norm(x)
    # z is just perpendicular to x and y
    z = np.cross(x, y)
    z /= np.linalg.norm(z)
    rotation = np.array([x, y, z]).T

    test = transform_residue(res, TranslationThenRotation(translation, rotation))
    assert np.all(test.get_atom(n.name).coord ** 2 < 0.001)
    assert np.all(test.get_atom(c1.name).coord[0:] < 0.001)
    assert -1.9 < test.get_atom(c1.name).coord[1] < -1.2, f"Unexpected C1' coordinates (x is bad): {test.get_atom(c1.name).coord}"
    assert np.all(test.get_atom(c2.name).coord[0] < -0.01)
    assert np.all(test.get_atom(c2.name).coord[1] > 0.01)

    return TranslationThenRotation(translation, rotation)

def fit_trans_rot_to_pairs(atom_pairs: Sequence[Tuple[Optional[Bio.PDB.Atom.Atom], Optional[Bio.PDB.Atom.Atom]]]) -> TranslationRotationTranslation:
    """
    Find a translation plus rotation minimizing the RMSD beweteen the atom pairs.
    """
    atom_pairs = [ (a, b) for a, b in atom_pairs if a is not None and b is not None ]
    if len(atom_pairs) < 3:
        raise ResideTransformError(f"Not enough pairs for fitting: {[ (a.id, b.id) for a, b in atom_pairs]}")
    coord_pairs = [ (a.coord, b.coord) for a, b in atom_pairs ]
    m1, m2 = zip(*coord_pairs)
    trans1 = -np.mean(m1, axis=0)
    m1 = np.array(m1) + trans1
    trans2 = np.mean(m2, axis=0)
    m2 = np.array(m2) - trans2

    rotation, rssd, _sensitivity_matrix = scipy.spatial.transform.Rotation.align_vectors(m1, m2, return_sensitivity=True) # type:ignore
    #                      return_sensitivity=True makes it throw on some unwanted special cases ^^^^^^^^^^^^^^^^^^

    return TranslationRotationTranslation(trans1, rotation.as_matrix(), trans2)

def atom_pair_rmsd(atom_pairs: Sequence[Tuple[Optional[Bio.PDB.Atom.Atom], Optional[Bio.PDB.Atom.Atom]]]) -> float:
    """
    Calculate RMSD between to sets of atoms, without rotating/translating them
    """
    atom_pairs = [ (a, b) for a, b in atom_pairs if a is not None and b is not None ]
    coord_pairs = [ (a.coord, b.coord) for a, b in atom_pairs ]
    return float(np.mean([ _linalg_norm(c1 - c2) for c1, c2 in coord_pairs ]))

@dataclass
class PairInformation:
    type: pair_defs.PairType
    res1: AltResidue
    res2: AltResidue
    position1: Optional[ResiduePosition]
    position2: Optional[ResiduePosition]
    rot_trans1: Optional[TranslationThenRotation]
    rot_trans2: Optional[TranslationThenRotation]

def calc_pair_information(pair_type: pair_defs.PairType, res1: AltResidue, res2: AltResidue) -> PairInformation:
    p1 = try_get_residue_posinfo(res1)
    p2 = try_get_residue_posinfo(res2)
    rot_trans1 = try_x(get_residue_posinfo_C1_N, res1)
    rot_trans2 = try_x(get_residue_posinfo_C1_N, res2)

    return PairInformation(pair_type, res1, res2, p1, p2, rot_trans1, rot_trans2)

def load_pair_information_by_idtuple(pair_type: Union[tuple[str, str], pair_defs.PairType], identifier: tuple[str, int, str, str, int, str, str, str, str, int, str, str]) -> PairInformation:
    pdbid, model, chain1, resname1, nr1, ins1, alt1, chain2, resname2, nr2, ins2, alt2 = identifier
    pair_type = pair_defs.PairType.from_tuple(pair_type)
    print(f"Loading ideal basepair:", *identifier)
    structure = pdb_utils.load_pdb(None, pdbid)

    res1 = get_residue(structure, None, model, chain1, resname1, nr1, ins1, alt1, None)
    res2 = get_residue(structure, None, model, chain2, resname2, nr2, ins2, alt2, None)
    # detach parent to lower memory usage (all structures would get replicated to all workers)
    res1 = res1.copy()
    res2 = res2.copy()
    assert res1.res.resname == resname1 and res2.res.resname == resname2, f"Residue names don't match: {res1.res.resname} {res2.res.resname} {resname1} {resname2} ({identifier})"
    return calc_pair_information(pair_type, res1, res2)

def load_ideal_pairs(pool: Union[Pool, MockPool], metadata_file: str) -> dict[pair_defs.PairType, PairInformation]:
    with open(metadata_file) as f:
        metadata = json.load(f)
    result = dict()
    
    for pair_meta in metadata:
        if not pair_meta.get("statistics", None):
            continue
        pair_type = pair_defs.PairType.from_tuple(pair_meta["pair_type"])
        if pair_type.n:
            continue
        largest_stat = max(pair_meta["statistics"], key=lambda x: x["count"])
        nicest_bp = largest_stat["nicest_bp"]
        result[pair_type] = pool.apply_async(load_pair_information_by_idtuple, args=[pair_type, nicest_bp])

    for k in result:
        result[k] = result[k].get()
    return result

class PairMetric: # "interface"
    columns=None
    def get_columns(self) -> Sequence[str]:
        if self.columns is None:
            raise NotImplementedError()
        return tuple(self.columns)

    def get_values(self, pair: PairInformation) -> Sequence[Optional[float]]:
        raise NotImplementedError()

    def select(self, columns: Sequence[str], rename_as: Optional[Sequence[str]]=None) -> 'PairMetric':
        metric = self
        class SelectPairMetric(PairMetric):
            def get_columns(self) -> Sequence[str]:
                return tuple(rename_as) if rename_as else tuple(columns)
            def get_values(self, pair: PairInformation) -> Sequence[Optional[float]]:
                selection = [ i for i, col in enumerate(metric.get_columns()) if col in columns ]
                vals = list(metric.get_values(pair))
                return tuple([vals[i] for i in selection])
        return SelectPairMetric()
    
    def apply_with(self, metric2: 'PairMetric', fn: Callable[[float, float], float], rename_as: Optional[Sequence[str]]=None) -> 'PairMetric':
        metric = self
        assert len(metric.get_columns()) == len(metric2.get_columns()), f"Cannot apply {fn.__name__} to metrics with different number of columns."

        class ApplyWithPairMetric(PairMetric):
            def get_columns(self) -> Sequence[str]:
                return tuple(rename_as) if rename_as else metric.get_columns()
            def get_values(self, pair: PairInformation) -> Sequence[Optional[float]]:
                return [ (None if a is None or b is None else fn(a, b))
                         for a, b in zip(metric.get_values(pair), metric2.get_values(pair))]

        return ApplyWithPairMetric()
    
    def apply(self, fn: Callable[[float], float], rename_as: Optional[Sequence[str]]=None) -> 'PairMetric':
        metric = self
        class ApplyPairMetric(PairMetric):
            def get_columns(self) -> Sequence[str]:
                return tuple(rename_as) if rename_as else metric.get_columns()
            def get_values(self, pair: PairInformation) -> Sequence[Optional[float]]:
                return [ (None if a is None else fn(a))
                         for a in metric.get_values(pair)]
        return ApplyPairMetric()

    def agg(self, fn: Callable[[list[Optional[float]]], float], rename_as: Optional[str]=None) -> 'PairMetric':
        metric = self
        class AggPairMetric(PairMetric):
            def get_columns(self) -> Sequence[str]:
                return (rename_as,) if rename_as else (metric.get_columns()[0],)
            def get_values(self, pair: PairInformation) -> Sequence[Optional[float]]:
                vals = list(metric.get_values(pair))
                return [fn(vals)]
        return AggPairMetric()

    @staticmethod
    def concat(metrics: list['PairMetric']) -> 'PairMetric':
        class ConcatPairMetric(PairMetric):
            def get_columns(self) -> Sequence[str]:
                return list(itertools.chain(*[m.get_columns() for m in metrics]))
            def get_values(self, pair: PairInformation) -> Sequence[Optional[float]]:
                return list(itertools.chain(*[m.get_values(pair) for m in metrics]))
        return ConcatPairMetric()

def Rx(θ): return np.array([ [1, 0, 0], [0, np.cos(θ), -np.sin(θ)], [0, np.sin(θ), np.cos(θ) ]])
def Rz(θ): return np.array([[ np.cos(θ), -np.sin(θ), 0],[ np.sin(θ), np.cos(θ), 0],[ 0, 0, 1] ])
def Ry(θ): return np.array([[ np.cos(θ), 0, np.sin(θ)],[ 0, 1, 0],[ -np.sin(θ), 0, np.cos(θ)]])

def get_C1_N_yaw_pitch_roll(rot1: np.ndarray, rot2: np.ndarray) -> Tuple[float, float, float]:
    matrix = rot1.T @ rot2
    # https://stackoverflow.com/a/37558238
    # yaw = math.degrees(math.atan2(matrix[1, 0], matrix[0, 0]))
    # pitch = math.degrees(math.atan2(-matrix[2, 0], math.sqrt(matrix[2, 1]**2 + matrix[2, 2]**2)))
    # roll = math.degrees(math.atan2(matrix[2, 1], matrix[2, 2]))

    # this function probably the intuitive representation, as it corresponds to:
    # In [57]: Rotation.from_matrix(Rz(0.33) @ Ry(0.44) @ Rx(0.55)).as_euler("ZYX")
    # Out[57]: array([0.33, 0.44, 0.55])
    yaw, pitch, roll = np.degrees(Rotation.from_matrix(matrix).as_euler("ZYX"))

    assert np.allclose(matrix, Rz(math.radians(yaw)) @ Ry(math.radians(pitch)) @ Rx(math.radians(roll)), atol=1e-5)
    return yaw, pitch, roll
def get_C1_N_euler_angles(rot1: np.ndarray, rot2: np.ndarray) -> Tuple[float, float, float]:
    matrix = rot1.T @ rot2
    # https://en.wikipedia.org/wiki/Rotation_formalisms_in_three_dimensions#Rotation_matrix_%E2%86%92_Euler_angles_(z-x-z_extrinsic)
    # eu1 = math.degrees(math.atan2(matrix[2, 0], matrix[2, 1]))
    # eu2 = math.degrees(math.acos(matrix[2, 2]))
    # eu3 = -math.degrees(math.atan2(matrix[0, 2], matrix[1, 2]))
    eu1, eu2, eu3 = np.degrees(Rotation.from_matrix(matrix).as_euler("zxz"))
    return eu1, eu2, eu3

T = TypeVar('T')
def try_x(f: Callable[[AltResidue], T], res: AltResidue, warning = None) -> Optional[T]:
    try:
        return f(res)
    except ResideTransformError as e:
        if warning is not None:
            warning(str(e))
        else:
            print(f"WARNING(residue = {res.res.full_id} {res.alt}): {e}")
        return None

def try_get_residue_posinfo(res: AltResidue, warning = None) -> Optional[ResiduePosition]:
    return try_x(get_residue_posinfo, res, warning)

TResidue = TypeVar('TResidue', bound=Union[AltResidue, Bio.PDB.Residue.Residue])
def transform_residue(r: TResidue, t: Union[TranslationThenRotation, TranslationRotationTranslation, ResiduePosition]) -> TResidue:
    if isinstance(r, AltResidue):
        return AltResidue(transform_residue(r.res, t), r.alt) # type: ignore
    r = r.copy()
    if isinstance(t, TranslationThenRotation):
        r.transform(np.eye(3), t.translation)
        r.transform(t.rotation, np.zeros(3))
    elif isinstance(t, ResiduePosition):
        r.transform(np.eye(3), t.origin)
        r.transform(t.rotation, np.zeros(3))
    elif isinstance(t, TranslationRotationTranslation):
        r.transform(np.eye(3), t.translation1)
        r.transform(t.rotation, np.zeros(3))
        r.transform(np.eye(3), t.translation2)
    return r

def to_float32_df(df: pl.DataFrame) -> pl.DataFrame:
    """Convert float64 columns to float32"""
    return df.with_columns([
        pl.col(c).cast(pl.Float32).alias(c) for c in df.columns if df[c].dtype == pl.Float64
    ])

def _linalg_norm(v) -> float:
    return float(np.linalg.norm(v))

def _linalg_normalize(v) -> np.ndarray:
    return v / np.linalg.norm(v)

def maybe_agg(agg, array):
    """
    Apply a function to a list if it's not empty
    """
    if len(array) == 0:
        return None
    return agg(array)

def abs_min(array):
    """
    Finds the element closest to zero
    """
    ix = np.argmin(np.abs(array))
    return array[ix]

def get_edge_atoms(res: AltResidue, edge_name: str | Literal['W', 'H', 'S']) -> list[Bio.PDB.Atom.Atom]:
    """
    Returns the residue atoms forming a given edge
    """
    edge = pair_defs.base_edges.get(resname_map.get(res.resname, res.resname), {}).get(edge_name, [])
    return [ a for aname in edge if (a:=res.get_atom(aname, None)) is not None ]

def get_edge_edges(edge: list[Bio.PDB.Atom.Atom]) -> Optional[tuple[Bio.PDB.Atom.Atom, Bio.PDB.Atom.Atom]]:
    """
    Returns the first and last atom of an edge
    """
    if not edge or len(edge) < 2:
        return None
    polar_atoms = [ a for a in edge if a.element not in ("H", "C") ]
    if len(polar_atoms) < 2:
        return (edge[0], edge[-1])
    return (polar_atoms[0], polar_atoms[-1])

class CoplanarityEdgeMetrics(PairMetric):
    columns = [ "coplanarity_angle", "coplanarity_shift1", "coplanarity_shift2", "coplanarity_edge_angle1", "coplanarity_edge_angle2" ]
    def get_values(self, pair: PairInformation) -> Sequence[Optional[float]]:
        p1, p2 = pair.position1, pair.position2
        if p1 is None or p2 is None:
            return [ None ] * len(self.columns)

        edge1 = get_edge_atoms(pair.res1, pair.type.edge1_name)
        edge1_edges = get_edge_edges(edge1)
        edge2 = get_edge_atoms(pair.res2, pair.type.edge2_name)
        edge2_edges = get_edge_edges(edge2)

        return [
            math.degrees(np.arccos(np.dot(p1.plane_normal_vector, p2.plane_normal_vector))), # coplanarity_angle
            maybe_agg(abs_min, [p2.plane_get_deviation(a.coord) for a in edge1]), # coplanarity_shift1
            maybe_agg(abs_min, [p1.plane_get_deviation(a.coord) for a in edge2]), # coplanarity_shift2
            p2.get_relative_line_rotation(*(a.coord for a in edge1_edges)) if edge1_edges is not None else None, # coplanarity_edge_angle1
            p1.get_relative_line_rotation(*(a.coord for a in edge2_edges)) if edge2_edges is not None else None, # coplanarity_edge_angle2
        ]
    
class IsostericityCoreMetrics(PairMetric):
    """
    C1'-C1' distance and the angle of the two C1'-N bonds, as the two numbers reported in the LSW2002 paper
    """
    columns = [ "C1_C1_distance", "C1_C1_total_angle" ]
    def get_values(self, pair: PairInformation) -> Sequence[Optional[float]]:
        pos1, pos2 = pair.position1, pair.position2
        if pos1 is None or pos2 is None or pos1.c1_pos is None or pos2.c1_pos is None:
            return [ None ] * len(self.columns)

        distance = _linalg_norm(pos1.c1_pos - pos2.c1_pos)
        if pos1.n_pos is None or pos2.n_pos is None:
            total_angle = None
        else:
            total_angle = math.degrees(np.arccos(np.dot(_linalg_normalize(pos1.n_pos - pos1.c1_pos), _linalg_normalize(pos2.c1_pos - pos2.n_pos))))
        return [ distance, total_angle ]
    
class EulerAngleMetrics(PairMetric):
    def __init__(self, col_names: Sequence[str], function: Callable[[np.ndarray, np.ndarray], tuple[float, float, float]], reference: Literal[1, 2]) -> None:
        self.columns = [ "C1_C1_" + c for c in col_names ]
        self.function = function
        self.reference = reference

    def get_values(self, pair: PairInformation) -> Sequence[Optional[float]]:
        rot_trans1, rot_trans2 = pair.rot_trans1, pair.rot_trans2
        if rot_trans1 is None or rot_trans2 is None:
            return [ None ] * len(self.columns)
        flip = np.diag([-1, -1, 1])
        if self.reference == 2:
            return self.function(rot_trans2.rotation, rot_trans1.rotation @ flip)
        else:
            return self.function(rot_trans1.rotation, rot_trans2.rotation @ flip)

class TranslationMetrics(PairMetric):
    """
    x, y, z coordinates in the other residue C1' - N1/N9 coordinate system
    """
    def __init__(self, reference: Literal[1, 2]):
        self.reference = reference
        self.columns = np.char.add(["x", "y", "z"], str(reference))
    def get_values(self, pair: PairInformation) -> Sequence[float | None]:
        ref, other = (pair.rot_trans1, pair.rot_trans2) if self.reference == 1 else (pair.rot_trans1, pair.rot_trans2)
        if ref is None or other is None:
            return [ None ] * 3
        return list(ref.rotation @ (other.translation - ref.translation))

class StandardMetrics:
    Coplanarity = CoplanarityEdgeMetrics()
    Isostericity = IsostericityCoreMetrics()
    EulerClassic = EulerAngleMetrics(["euler_phi", "euler_theta", "euler_psi"], get_C1_N_euler_angles, 1)
    YawPitchRoll = EulerAngleMetrics(["yaw1", "pitch1", "roll1"], get_C1_N_yaw_pitch_roll, 1)
    YawPitchRoll2 = EulerAngleMetrics(["yaw2", "pitch2", "roll2"], get_C1_N_yaw_pitch_roll, 2)
    Translation1 = TranslationMetrics(1)
    Translation2 = TranslationMetrics(2)

class RMSDToIdealMetric(PairMetric):
    """
    Various RMSD distances between the singular "ideal" basepair the observed one
    """
    def __init__(self, name, ideal: dict[pair_defs.PairType, PairInformation],
        fit_on: Literal['all', 'left_C1N', 'left', 'right', 'both_edges'],
        calculate: Literal['all', 'right_C1N', 'both_C1N', 'left_edge', 'right_edge', 'both_edges'],
    ) -> None:
        self.ideal = ideal
        self.columns = [ f"rmsd_{name}" ]
        self.fit_on = fit_on
        self.calculate = calculate

    def get_c1n_atoms(self, res: AltResidue) -> list[Bio.PDB.Atom.Atom]:
        n = res.get_atom("N9", None) or res.get_atom("N1", None)
        c1 = res.get_atom("C1", None)
        if n.name == "N9":
            # purine
            return [ c1, n, res.get_atom("C8", None), res.get_atom("C4", None) ] 
        else:
            # pyrimidine
            return [ c1, n, res.get_atom("C6", None), res.get_atom("C2", None) ]

    def atom_pairs(self, l, pair: PairInformation, ideal: PairInformation) -> list[tuple[Optional[Bio.PDB.Atom.Atom], Optional[Bio.PDB.Atom.Atom]]]:
        if l == 'all':
            return [ *self.atom_pairs('left', pair, ideal), *self.atom_pairs('right', pair, ideal) ]
        elif l == 'left':
            return [ (pair.res1.get_atom(a.name, None), ideal.res1.get_atom(a.name, None)) for a in get_base_atoms(pair.res1)[1] ]
        elif l == 'right':
            return [ (pair.res2.get_atom(a.name, None), ideal.res2.get_atom(a.name, None)) for a in get_base_atoms(pair.res2)[1] ]
        elif l == 'left_C1N':
            return list(zip(self.get_c1n_atoms(pair.res1), self.get_c1n_atoms(ideal.res1)))
        elif l == 'right_C1N':
            return list(zip(self.get_c1n_atoms(pair.res2), self.get_c1n_atoms(ideal.res2)))
        elif l == 'both_C1N':
            return [ *self.atom_pairs('left_C1N', pair, ideal), *self.atom_pairs('right_C1N', pair, ideal) ]
        elif l == 'left_edge' or l == 'right_edge':
            edge = pair.type.edge1_name if l == 'left_edge' else pair.type.edge2_name
            p_edge = get_edge_atoms(pair.res1 if l == 'left_edge' else pair.res2, edge)
            i_edge = get_edge_atoms(ideal.res1 if l == 'left_edge' else ideal.res2, edge)
            both_atoms = { a.name for a in p_edge } & { a.name for a in i_edge }
            return list(zip(
                [ a for a in p_edge if a.name in both_atoms ],
                [ a for a in i_edge if a.name in both_atoms ]
            ))
        elif l == 'both_edges':
            return [ *self.atom_pairs('left_edge', pair, ideal), *self.atom_pairs('right_edge', pair, ideal) ]
        else:
            raise NotImplementedError(f"atom_pairs={l} is not implemented")

    def get_values(self, pair: PairInformation) -> Sequence[Optional[float]]:
        ideal = self.ideal.get(pair.type.without_n(), None)
        if ideal is None or pair.rot_trans1 is None or pair.rot_trans2 is None:
            return [ None ]
        assert ideal.rot_trans1 is not None

        if self.fit_on == 'left_C1N':
            fit = TranslationRotationTranslation(pair.rot_trans1.translation, pair.rot_trans1.rotation @ ideal.rot_trans1.rotation.T, -ideal.rot_trans1.translation)
        else:
            atom_pairs = self.atom_pairs(self.fit_on, pair, ideal)
            assert all(a.name == b.name for (a, b) in atom_pairs if a and b)
            fit = fit_trans_rot_to_pairs(atom_pairs)
        # elif self.fit_on == 'left':
        #     atom_pairs = [ (pair.res1.get_atom(a, None), ideal.res1.get_atom(a, None)) for a in get_base_atoms(pair.res1)[1] ]
        #     fit = fit_trans_rot_to_pairs(atom_pairs)
        # else:
        #     raise NotImplementedError(f"fit_on={self.fit_on}")

        transformed1 = transform_residue(pair.res1, fit)
        transformed2 = transform_residue(pair.res2, fit)

        calc_atom_pairs = self.atom_pairs(self.calculate, dataclasses.replace(pair, res1=transformed1, res2=transformed2), ideal)
        assert all(a.name == b.name for (a, b) in calc_atom_pairs if a and b)
        rmsd = atom_pair_rmsd(calc_atom_pairs)
        return [ rmsd ]


def get_angle(atom1: Bio.PDB.Atom.Atom, atom2: Bio.PDB.Atom.Atom, atom3: Bio.PDB.Atom.Atom) -> Optional[float]:
    """3 point angle (in degrees)"""
    if atom1 is None or atom2 is None or atom3 is None:
        return None
    v1 = atom1.coord - atom2.coord
    v2 = atom3.coord - atom2.coord
    return math.degrees(np.arccos(np.dot(v1, v2) / (_linalg_norm(v1) * _linalg_norm(v2))))

def get_distance(atom1: Bio.PDB.Atom.Atom, atom2: Bio.PDB.Atom.Atom) -> float:
    return _linalg_norm(atom1.coord - atom2.coord)


@dataclass
class HBondStats:
    length: Optional[float]
    donor_angle: Optional[float]
    acceptor_angle: Optional[float]
    # "out-of-plane angle" - angle between the line of the bond and the first base plane
    OOPA1: Optional[float] = None
    OOPA2: Optional[float] = None

def hbond_stats(
    atom0: Bio.PDB.Atom.Atom, atom1: Bio.PDB.Atom.Atom, atom2: Bio.PDB.Atom.Atom, atom3: Bio.PDB.Atom.Atom,
    residue1_plane: Optional[ResiduePosition],
    residue2_plane: Optional[ResiduePosition]
) -> HBondStats:
    """Gets stats for a given hydrogen bond in a given basepair."""
    if atom1 is None or atom2 is None:
        return HBondStats(None, None, None)

    assert _linalg_norm(atom0.coord - atom1.coord) <= 1.6, f"atoms 0,1 not bonded: {atom0.full_id} {atom1.full_id} ({np.linalg.norm(atom0.coord - atom1.coord)} > 1.6, {atom0.coord} {atom1.coord})"
    assert _linalg_norm(atom2.coord - atom3.coord) <= 1.6, f"atoms 2,3 not bonded: {atom2.full_id} {atom3.full_id} ({np.linalg.norm(atom0.coord - atom1.coord)} > 1.6, {atom0.coord} {atom1.coord})"
    assert _linalg_norm(atom1.coord - atom2.coord) > 1.55, f"atoms too close for h-bond: {atom1.full_id} {atom2.full_id} ({np.linalg.norm(atom1.coord - atom2.coord)} < 2, {atom1.coord} {atom2.coord})"
    assert _linalg_norm(atom1.coord - atom2.coord) < 20, f"atoms too far for h-bond: ({np.linalg.norm(atom1.coord - atom2.coord)} > 6, {atom1.coord} {atom2.coord})"

    result = HBondStats(
        length=get_distance(atom1, atom2),
        donor_angle=get_angle(atom0, atom1, atom2),
        acceptor_angle=get_angle(atom3, atom2, atom1)
    )

    if residue1_plane is not None:
        result.OOPA1 = residue1_plane.get_relative_line_rotation(atom1.coord, atom2.coord)
    if residue2_plane is not None:
        result.OOPA2 = residue2_plane.get_relative_line_rotation(atom2.coord, atom1.coord)

    return result

resname_map = {
    'DT': 'U',
    'DC': 'C',
    'DA': 'A',
    'DG': 'G',
    'DU': 'U',
    'T': 'U',
}

def pair_hbonds_stats(pair: PairInformation) -> Optional[List[Optional[HBondStats]]]:
    """Returns stats for all h-bonds in the given basepair (assuming the `pair.type`)"""
    r1, r2 = pair.res1, pair.res2
    rp1, rp2 = pair.position1, pair.position2
    # if r1.resname > r2.resname:
    #     r1, r2 = r2, r1
    hbonds = pdef.get_hbonds(pair.type, throw=False)
    if len(hbonds) == 0:
        return None
    # print(pair_type, pair_name, hbonds)
    def get_atom(n):
        if n[0] == 'A':
            return r1.get_atom(n[1:], None)
        elif n[0] == 'B':
            return r2.get_atom(n[1:], None)
        else:
            raise ValueError(f"Invalid atom name: {n}")

    bonds_atoms = [
        tuple(get_atom(a) for a in atoms)
        for atoms in hbonds
    ]

    return [
        hbond_stats(atoms[0], atoms[1], atoms[2], atoms[3], rp1, rp2)
        for atoms, hb in zip(bonds_atoms, hbonds)
    ]

def get_atom_df(model: Bio.PDB.Model.Model) -> dict[str, np.ndarray]:
    atom_count = sum(1 for atom in model.get_atoms())
    atom_res = np.zeros(atom_count, dtype=np.int32)
    atom_resname = np.zeros(atom_count, dtype="S4")
    atom_chain = np.zeros(atom_count, dtype="S4")
    atom_coord = np.zeros((atom_count, 3), dtype=np.float32)
    atom_element = np.zeros(atom_count, dtype="S4")
    atom_name = np.zeros(atom_count, dtype="S4")

    atom_idx = 0
    for chain in model.get_chains():
        for residue in chain.get_residues():
            for atom in residue.get_atoms():

                atom_res[atom_idx] = residue.id[1]
                atom_resname[atom_idx] = residue.resname
                atom_chain[atom_idx] = chain.id
                atom_coord[atom_idx, :] = atom.coord
                atom_element[atom_idx] = atom.element or ""
                atom_name[atom_idx] = (atom.name or '').encode('utf-8')

                atom_idx += 1

    return {
        "chain": atom_chain,
        "res": atom_res,
        "resname": atom_resname,
        "coord": atom_coord,
        "element": atom_element,
        "name": atom_name,
    }

def to_csv_row(df: pl.DataFrame, i: int = 0, max_len = None) -> str:
    row = next(df[i, :].iter_rows())
    if max_len is not None:
        row = row[:max_len]
    return ','.join(str(x) for x in row)

T = TypeVar("T")
def lazy(create: Callable[..., T], *args) -> Callable[[], T]:
    cache: T
    has_value = False
    def core():
        nonlocal cache, has_value
        if not has_value:
            cache = create(*args)
            has_value = True
        return cache
    return core

def get_residue(structure: Bio.PDB.Structure.Structure, strdata: Optional[Callable[[], pdb_utils.StructureData]], model, chain, resname, nr, ins, alt, symop) -> AltResidue:
    """Finds a residue in a given structure, and optionally performs the given symmetry operation on it."""
    nr = int(nr)
    ins = ins or ' '
    alt = alt or ''
    ch: Bio.PDB.Chain.Chain = structure[model-1][chain]
    found: Bio.PDB.Residue.Residue = ch.child_dict.get((' ', nr, ins), None)
    if found is None:
        found = next((r for r in ch.get_residues() if r.id[1] == nr and r.id[2] == ins), None)
    if found is None:
        raise KeyError(f"{chain}.{nr}{ins}")
    
    if symop:
        assert strdata is not None
        sym_transform = next((a for a in strdata().assembly if a.pdbname == symop or f"ASM_{a.id}" == symop), None)
        if sym_transform is not None:
            found = found.copy()
            found.transform(sym_transform.rotation, sym_transform.translation)
        else:
            raise KeyError(f"Symmetry operation {symop} is not defined {structure.id}")
        
    if isinstance(found, Bio.PDB.Residue.DisorderedResidue):
        disordered_residues: list[Bio.PDB.Residue.Residue] = list(found.child_dict.values())
        disordered_residues.sort(key=lambda x: (resname == x.resname, len([ 0 for a in x.get_atoms() if alt in a.disordered_get_id_list() ])))
        found = disordered_residues[-1]

    return AltResidue(found, alt)


def get_stats_for_csv(df_: pl.DataFrame, structure: Bio.PDB.Structure.Structure, strdata: Callable[[], pdb_utils.StructureData], metrics: list[PairMetric]) -> Iterator[tuple[int, list[Optional[HBondStats]], list[Optional[float]]]]:
    """
    Iterates over all basepairs in the given DataFrame, calculating the stats for each one. Works on a single PDB structure.
    """
    df = df_.with_row_count()
    pair_type_col = df["family" if "family" in df.columns else "type"]
    metric_columns = [ tuple(m.get_columns()) for m in metrics ]
    for (pair_family, i, pdbid, model, chain1, resname1, nr1, ins1, alt1, symop1, chain2, resname2, nr2, ins2, alt2, symop2) in zip(pair_type_col, df["row_nr"], df["pdbid"], df["model"], df["chain1"], df['res1'], df["nr1"], df["ins1"], df["alt1"], df["symmetry_operation1"], df["chain2"], df['res2'], df["nr2"], df["ins2"], df["alt2"], df["symmetry_operation2"]):
        try:
            assert structure.id.lower() == pdbid.lower(), f"pdbid mismatch: {structure.id} != {pdbid}"
            ins1 = (ins1 or '').strip()
            ins2 = (ins2 or '').strip()

            try:
                r1 = get_residue(structure, strdata, model, chain1, resname1, nr1, ins1, alt1, symop1)
            except KeyError as keyerror:
                print(f"Could not find residue1 {pdbid}.{model}.{chain1}.{str(nr1)+ins1} ({keyerror=})")
                continue
            try:
                r2 = get_residue(structure, strdata, model, chain2, resname2, nr2, ins2, alt2, symop2)
            except KeyError as keyerror:
                print(f"Could not find residue2 {pdbid}.{model}.{chain2}.{str(nr2)+ins2} in {pdbid} ({keyerror=})")
                continue
            pair_type = pair_defs.PairType.create(pair_family, r1.resname, r2.resname, name_map=resname_map)
            pair_info = calc_pair_information(pair_type, r1, r2)
            hbonds = pair_hbonds_stats(pair_info)
            metric_values: list[Optional[float]] = []
            for mix, m in enumerate(metrics):
                v = m.get_values(pair_info)
                assert len(v) == len(metric_columns[mix]), f"Invalid number of values from metric {m} (got {len(v)}, expected {metric_columns[mix]})"
                metric_values.extend(v)
            if hbonds is not None:
                yield i, hbonds, metric_values
        except AssertionError as e:
            print(f"Assertion error {e} on row:\n{to_csv_row(df_, i, 18)}")
            continue
        except Exception as keyerror:
            print(f"Error on row:\n{to_csv_row(df_, i, 18)}")
            import traceback
            print(keyerror)
            print(traceback.format_exc())
            print("Continuing...")
            continue

def remove_duplicate_pairs_phase1(df: pl.DataFrame):
    """
    Removes duplicate symmetric pairs - i.e. cHS-AG and cHS-GA.
    Keeps only the preferred family ordering (WS > SW), and only the preferred base ordering (GC > CG)
    If the basepair family has a lowercase letter, it will always be second (preserving the FR3D ordering).
    """
    # original = df
    original_len = len(df)
    df = df.with_row_index("_tmp_row_nr")
    family_col = pl.col("family" if "family" in df.columns else "type")
    df = df.with_columns(
        family_col.str.replace("^(n?[ct])([WHSBwhsb])([WHSBwhsb])([a-z]*)$", "$1$3$2$4").alias("_tmp_family_reversed")
    )
    def core(df: pl.DataFrame, condition):
        duplicated = df.filter(condition)\
            .join(df,
                left_on=[family_col, 'pdbid', 'model', 'chain1', 'res1', 'nr1', 'ins1', 'alt1', 'chain2', 'nr2', 'ins2', 'alt2'],
                right_on=['_tmp_family_reversed', 'pdbid', 'model', 'chain2', 'res2', 'nr2', 'ins2', 'alt2', 'chain1', 'nr1', 'ins1', 'alt1'],
                join_nulls=True, how="inner"
            )
        duplicated_f = duplicated.select("_tmp_row_nr", "_tmp_row_nr_right")
        duplicated_set = set(duplicated_f["_tmp_row_nr"])
        check_not_removed = set(duplicated_f["_tmp_row_nr_right"])
        if duplicated_set.intersection(check_not_removed):
            # print(f"{len(duplicated_set)} entries would be removed twice! {list(duplicated_set.intersection(check_not_removed))}")
            # print(df.filter(pl.col("_tmp_row_nr").is_in(duplicated_set.intersection(check_not_removed))).head(10))

            # for x in list(duplicated_set.intersection(check_not_removed))[:10]:
            #     print(x)
            #     print(duplicated.filter(pl.col("_tmp_row_nr") == x).select("family", pl.col("res1") + pl.col("res2"), "pdbid", "chain1", "chain2", pl.col("alt1") + pl.col("nr1").cast(pl.Utf8) + pl.col("ins1"), pl.col("alt2") + pl.col("nr2").cast(pl.Utf8) + pl.col("ins2")))
            #     print(duplicated.filter(pl.col("_tmp_row_nr_right") == x).select("family", pl.col("res1") + pl.col("res2"), "pdbid", "chain1", "chain2", pl.col("alt1") + pl.col("nr1").cast(pl.Utf8) + pl.col("ins1") , pl.col("alt2") + pl.col("nr2").cast(pl.Utf8) + pl.col("ins2")))
            # assert False

            duplicated_set.difference_update(check_not_removed)

        df = df.filter(pl.col("_tmp_row_nr").is_in(duplicated_set).not_())
        return df

    df = core(df, family_col.str.contains("^(?i)n?(c|t)(SH|HW|SW|BW)a?$"))
    # print(f"Removed duplicates - family orientation {original_len} -> {len(df)}")

    df = core(df, family_col.str.contains("^n?[ct][whsb][WHSB]a?$"))
    # print(f"Removed duplicates - FR3D small letter is second {original_len} -> {len(df)}")

    base_ordering = { 'A': 1, 'G': 2, 'C': 3, 'U': 4, 'T': 4 }
    df = core(df,
                family_col.str.to_lowercase().is_in(["tss", "css"]).not_() &
                (pl.col("res1").replace(resname_map).replace(base_ordering, default=0) > pl.col("res2").replace(resname_map).replace(base_ordering, default=0)))
    # print(f"Removed duplicates - base ordering {original_len} -> {len(df)}")

    assert len(df) >= original_len / 2

    return df.drop([ "_tmp_row_nr", "_tmp_family_reversed" ])

def remove_duplicate_pairs(df: pl.DataFrame):
    """
    Removes duplicate symmetric pairs - i.e. cHS-AG and cHS-GA.
    If the order is ambiguous according to `remove_duplicate_pairs_phase1`, keeps the one with shorter h-bonds, or with lower left residue number
    """
    df = remove_duplicate_pairs_phase1(df)

    def pair_id(type, res1, res2, chain1, nr1, ins1, alt1, chain2, nr2, ins2, alt2):
        pair = [ (chain1, nr1, ins1, alt1), (chain2, nr2, ins2, alt2) ]
        pair.sort()
        bases = resname_map.get(res1, res1) + "-" + resname_map.get(res2, res2)
        return type + "-" + bases + "|".join(tuple(str(x) for x in (pair[0] + pair[1])))
    family_col = "family" if "family" in df.columns else "type"
    pair_ids = [ pair_id(type, res1, res2, chain1, nr1, ins1, alt1, chain2, nr2, ins2, alt2) for type, res1, res2, chain1, nr1, ins1, alt1, chain2, nr2, ins2, alt2 in zip(df[family_col], df["res1"], df["res2"], df["chain1"], df["nr1"], df["ins1"], df["alt1"], df["chain2"], df["nr2"], df["ins2"], df["alt2"]) ]
    df = df.with_columns(
        pl.Series(pair_ids, dtype=pl.Utf8).alias("_tmp_pair_id")
    ).with_row_count("_tmp_row_nr")
    score = pl.lit(0, dtype=pl.Float64)
    null_count = pl.lit(0, dtype=pl.Int64)
    for col in df.columns:
        if col.startswith("hb_") and col.endswith("_length"):
            score += pl.col(col).fill_null(0)
            null_count += pl.col(col).is_null().cast(pl.Int64)
        if col.startswith("dssr_"):
            score += pl.col(col).is_null().cast(pl.Float64) * 0.03

    sym_pt = [ x for x in pair_defs.defined_pair_types() if pair_defs.has_symmetrical_definition(x) ]
    is_symmetrical_definition = (pl.col("family") + "-" + pl.col("res1").replace(resname_map) + "-" + pl.col("res2").replace(resname_map)).str.to_lowercase().is_in(set(sym_pt.__str__().lower()))

    score = pl.when(is_symmetrical_definition).then(pl.lit(0.0)).otherwise(score)

    df = df.sort([null_count, score, "chain1", "nr1", "ins1", "alt1", "chain2", "nr2", "ins2", "alt2"], descending=False)
    df = df.unique([family_col, "pdbid", "model", "_tmp_pair_id"], keep="first", maintain_order=True)
    df = df.sort("_tmp_row_nr")
    df = df.drop(["_tmp_pair_id", "_tmp_row_nr"])
    return df

def get_max_bond_count(df: pl.DataFrame):
    """Gets the maximum number of hydrogen bonds in the given DataFrame.  """
    pair_types = list(set(pair_defs.PairType.create(type, res1, res2, name_map=resname_map) for type, res1, res2 in df[["family", "res1", "res2"]].unique().iter_rows()))
    bond_count = max(len(pair_defs.get_hbonds(p, throw=False)) for p in pair_types)
    # print("Analyzing pair types:", pair_types, "with bond count =", bond_count)
    print(f"Analyzing {len(pair_types)} pair  typeswith bond count =", bond_count)
    return max(3, bond_count)

def make_backbone_columns(df: pl.DataFrame, structure: Optional[Bio.PDB.Structure.Structure]) -> list[pl.Series]:
    """Calculates is_dinucleotide and is_parallel columns"""
    is_dinucleotide = pl.Series("is_dinucleotide", [ False ] * len(df), dtype=pl.Boolean)
    is_parallel = pl.Series("is_parallel", [ None ] * len(df), dtype=pl.Boolean)

    if structure is None:
        return [ is_dinucleotide, is_parallel ]

    chain_index = {
        (model.serial_num, chain.id): {
            (res.id[1], res.id[2]): i
            for i, res in enumerate(chain.get_residues())
            if res.resname != 'HOH'
        }
        for model in structure.get_models()
        for chain in model.get_chains()
    }

    pair_set = set(
        (model, chain1, chain_index[(model, chain1)][(nr1, ins1 or ' ')], chain2, chain_index[(model, chain2)][(nr2, ins2 or ' ')])
        for model, chain1, nr1, ins1, chain2, nr2, ins2 in zip(df["model"], df["chain1"], df["nr1"], df["ins1"], df["chain2"], df["nr2"], df["ins2"])
        if  (nr1, ins1 or ' ') in chain_index.get((model, chain1), []) and (nr2, ins2 or ' ') in chain_index.get((model, chain2), [])
    )

    for i, (pdbid, model, chain1, nr1, ins1, alt1, chain2, nr2, ins2, alt2) in enumerate(zip(df["pdbid"], df["model"], df["chain1"], df["nr1"], df["ins1"], df["alt1"], df["chain2"], df["nr2"], df["ins2"], df["alt2"])):

        res_ix1 = chain_index.get((model, chain1), {}).get((nr1, ins1 or ' '), None)
        res_ix2 = chain_index.get((model, chain2), {}).get((nr2, ins2 or ' '), None)
        if res_ix1 is None or res_ix2 is None:
            continue
        is_dinucleotide[i] = abs(res_ix2 - res_ix1) == 1

        par = (model, chain1, res_ix1 + 1, chain2, res_ix2 + 1) in pair_set or \
           (model, chain1, res_ix1 - 1, chain2, res_ix2 - 1) in pair_set
        antipar = (model, chain1, res_ix1 + 1, chain2, res_ix2 - 1) in pair_set or \
           (model, chain1, res_ix1 - 1, chain2, res_ix2 + 1) in pair_set

        if antipar and not par:
            is_parallel[i] = False
        elif par and not antipar:
            is_parallel[i] = True

        # if antipar and par:
        #     print("Interesting pair:", pdbid, model, chain1, nr1, ins1, chain2, nr2, ins2, "both parallel and antiparallel")
    return [ is_dinucleotide, is_parallel ]



def make_stats_columns(pdbid: str, df: pl.DataFrame, add_metadata_columns: bool, max_bond_count: int, metrics: list[PairMetric]) -> Tuple[str, pl.DataFrame, np.ndarray]:
    """
    Adds the columns with basepair stats to the DataFrame.
    """
    bond_params = [ x.name for x in dataclasses.fields(HBondStats) ]
    valid = np.zeros(len(df), dtype=np.bool_)
    hb_columns: list[pl.Series] = [
        pl.Series(f"hb_{i}_{p}", [ None ] * len(df), dtype=pl.Float32)
        for i in range(max_bond_count)
        for p in bond_params ]
    metric_columns: list[pl.Series] = [
        pl.Series(name, [ None ] * len(df), dtype=pl.Float32)
        for m in metrics
        for name in m.get_columns()
    ]

    structure = None
    try:
        structure = pdb_utils.load_pdb(None, pdbid)
    except HTTPError as e:
        print(f"Could not load structure {pdbid} due to http error")
        print(e)
    except Exception as e:
        import traceback
        print(f"Could not load structure {pdbid} due to unknown error")
        print(e)
        print(traceback.format_exc())

    structure_data = lazy(lambda: pdb_utils.load_sym_data(None, pdbid))

    if structure is not None:
        for i, hbonds, metric_values in get_stats_for_csv(df, structure, structure_data, metrics):
            valid[i] = True
            for j, s in enumerate(hbonds):
                if s is None:
                    continue

                for param_i, param in enumerate(bond_params):
                    hb_columns[j * len(bond_params) + param_i][i] = getattr(s, param)

            for j, m in enumerate(metric_values):
                metric_columns[j][i] = m

    result_df = pl.DataFrame([*hb_columns, *metric_columns])
    result_df = to_float32_df(result_df)
    if add_metadata_columns:
        h = structure.header if structure is not None else dict()
        result_df = result_df.with_columns(
            pl.lit(h.get('deposition_date', None), dtype=pl.Utf8).alias("deposition_date"),
            pl.lit(h.get('name', None), dtype=pl.Utf8).alias("structure_name"),
            pl.lit(h.get('structure_method', None), dtype=pl.Utf8).alias("structure_method"),
            pl.lit(h.get('resolution', None), dtype=pl.Float32).alias("resolution")
        )

    result_df = result_df.with_columns(make_backbone_columns(df, structure))
    assert len(result_df) == len(df)
    return pdbid, result_df, valid

def df_hstack(columns: list[pl.DataFrame]) -> pl.DataFrame:
    return functools.reduce(pl.DataFrame.hstack, columns)

def load_input_pdbids(args) -> list[str]:
    r = set()
    for i in args.inputs:
        if i.endswith(".parquet"):
            r.update(pl.scan_parquet(i).select(pl.col("pdbid").str.to_lowercase()).unique().collect(streaming=True)["pdbid"].to_list())
        elif i.endswith(".csv"):
            r.update(pl.scan_parquet(i).select(pl.col("pdbid").str.to_lowercase()).unique().collect(streaming=True)["pdbid"].to_list())
        elif i.endswith("_basepair.txt") or i.endswith("_basepair_detail.txt"):
            r.add(os.path.basename(i).split("_")[0].lower())
        else:
            raise ValueError(f"Unknown input file type: {i}")
    return list(sorted(r))

def load_inputs(pool: Union[Pool, MockPool], args, pdbid_prefix='') -> pl.DataFrame:
    """Loads CSV/Parquet/FR3D files and returns a DataFrame with basepair identifiers"""
    inputs: list[str] = args.inputs
    if len(inputs) == 0:
        raise ValueError("No input files specified")
    if pdbid_prefix:
        pdbid_filter = pl.col('pdbid').str.to_lowercase().str.starts_with(pdbid_prefix.lower())
    else:
        pdbid_filter = pl.lit(True)
    
    if inputs[0].endswith(".parquet") or inputs[0].endswith('.csv'):
        print(f'Loading basepairing CSV files')
        df = scan_pair_csvs(args.inputs)
        df = df.filter(pdbid_filter)
        df = df.with_columns(
            alt1=pl.when((pl.col("alt1") == '\0') | (pl.col("alt1") == '?')).then(pl.lit('')).otherwise(pl.coalesce(pl.col("alt1"), pl.lit('')).str.strip_chars()),
            alt2=pl.when((pl.col("alt2") == '\0') | (pl.col("alt2") == '?')).then(pl.lit('')).otherwise(pl.coalesce(pl.col("alt2"), pl.lit('')).str.strip_chars()),
            ins1=pl.coalesce(pl.col("ins1"), pl.lit('')).str.strip_chars(),
            ins2=pl.coalesce(pl.col("ins2"), pl.lit('')).str.strip_chars(),
        )
        if "pdbsymstr" in df.columns:
            df = df.with_columns(
                symmetry_operation1=pl.lit(None, pl.Utf8),
                symmetry_operation2=pl.when(pl.col("pdbsymstr") == "1_555").then(pl.lit(None, pl.Utf8)).otherwise(pl.col("pdbsymstr"))
            )
        allowed_residues = ["A", "T", "G", "U", "C", "DA", "DT", "DG", "DC", "DU"]
        df = df.filter(pl.col("res1").is_in(allowed_residues)).filter(pl.col("res2").is_in(allowed_residues))
        df = df.sort('pdbid', 'model', 'nr1', 'nr2')
        # if "r11" in df.columns:
        df = df.collect(streaming=True)
    elif inputs[0].endswith("_basepair.txt") or inputs[0].endswith("_basepair_detail.txt"):
        print(f"Loading {len(inputs)} basepair files")
        import fr3d_parser
        df = fr3d_parser.read_fr3d_files_df(pool, args.inputs,
            filter=pdbid_filter,
            # filter=(pl.col("symmetry_operation1").is_null() & pl.col("symmetry_operation2").is_null())
        ).sort('pdbid', 'model', 'nr1', 'nr2')
    else:
        raise ValueError("Unknown input file type")
    if "type" not in df.columns and "family" not in df.columns:
        if not args.override_pair_family:
            raise ValueError("Input does not contain family column and --override-pair-type was not specified")
    if args.override_pair_family:
        df = df.select(
            pl.lit(args.override_pair_family).alias("type"),
            pl.lit(args.override_pair_family).alias("family"),
            pl.col("*"))
    return df

def validate_missing_columns(chunks: list[pl.DataFrame]):
    all_columns = set(col for chunk in chunks for col in chunk.columns)
    for c in chunks:
        if not all_columns.issubset(c.columns):
            print(f"Chunk with {set(c['pdbid'].to_numpy())} is missing columns: {all_columns.difference(set(c.columns))}")

def main_partition(pool: Union[Pool, MockPool], args, pdbid_partition_prefix='', ideal_basepairs: Optional[dict[pair_defs.PairType, PairInformation]] = None):
    df = load_inputs(pool, args, pdbid_prefix=pdbid_partition_prefix)
    if args.export_only:
        print("Exporting metadata only")
        return df
    raw_df_len = len(df)
    print(f"Loaded {raw_df_len} raw basepairs")
    if args.disable_cross_symmetry:
        df = df.filter(pl.col("symmetry_operation1").is_null() & pl.col("symmetry_operation2").is_null())
        print(f"Removed cross-symmetry basepairs -> len={len(df)}")
    else:
        if 'symmetry_operation1' in df.columns and 'symmetry_operation2' in df.columns:
            # at least one has to be in the primary asymmetric unit
            old_pdbids = set(df["pdbid"])
            df = df.filter(pl.col("symmetry_operation1").is_null() | pl.col("symmetry_operation2").is_null())
            if len(df) < raw_df_len:
                print(f"Removed purely cross-symmetry basepairs -> len={len(df)}, removed structs={old_pdbids - set(df['pdbid'])}")
            del old_pdbids
    if args.dedupe:
        df = remove_duplicate_pairs_phase1(df)
        print(f"Removed duplicates, phase1 -> len={len(df)}")
    max_bond_count = get_max_bond_count(df)

    groups = list(df.group_by(pl.col("pdbid")))

    processes = []

    pair_metrics = [
        StandardMetrics.Coplanarity,
        StandardMetrics.Isostericity,
        StandardMetrics.EulerClassic,
        StandardMetrics.YawPitchRoll,
        StandardMetrics.YawPitchRoll2,
        StandardMetrics.Translation1,
        StandardMetrics.Translation2,
    ]
    if ideal_basepairs:
        print(f"Loaded {len(ideal_basepairs)} ideal basepairs")
        pair_metrics.extend([
            RMSDToIdealMetric('C1N_frames1', ideal_basepairs, fit_on='left_C1N', calculate='right_C1N'),
            RMSDToIdealMetric('C1N_frames2', ideal_basepairs, fit_on='right_C1N', calculate='left_C1N'),
            RMSDToIdealMetric('edge1', ideal_basepairs, fit_on='left', calculate='right_edge'),
            RMSDToIdealMetric('edge2', ideal_basepairs, fit_on='right', calculate='left_edge'),
            RMSDToIdealMetric('edge_C1N_frame', ideal_basepairs, fit_on='both_edges', calculate='both_C1N'),
            RMSDToIdealMetric('all_base', ideal_basepairs, fit_on='all', calculate='all'),
        ])

    processes.append([ # process per PDB structure
        pool.apply_async(make_stats_columns, args=[pdbid, group, args.metadata, max_bond_count, pair_metrics])
        for (pdbid,), group in groups # type: ignore
    ])
    if args.dssr_binary is not None:
        import dssr_wrapper
        processes.append([
            pool.apply_async(dssr_wrapper.add_dssr_info, args=[pdbid, group, args.dssr_binary])
            for ((pdbid,), group) in groups # type: ignore
        ])

    def postfilter(df: pl.DataFrame):
        if args.postfilter_hb:
            df = df.filter(pl.any_horizontal(pl.col("^hb_\\d+_length$") < args.postfilter_hb))
        if args.postfilter_shift:
            df = df.filter((pl.col("coplanarity_shift1").abs() < args.postfilter_shift) & (pl.col("coplanarity_shift2").abs() < args.postfilter_shift))
        if args.dedupe:
            df = remove_duplicate_pairs(df)
        return df


    result_chunks = []
    for ((_pdbid,), group), ps in zip(groups, zip(*processes)):
        # each group of processes returns a tuple
        # 0. PDBID of the chunk
        # 1. a set of columns to be added to the original DataFrame
        # 2. a boolean array indicating which rows are valid or should be dropped
        ps = [ p.get() for p in ps ]
        pdbids = [ p[0] for p in ps ]
        assert pdbids == [ _pdbid ] * len(pdbids), f"pdbid mismatch: {_pdbid} x {pdbids}"
        added_columns, valid = [ x[1] for x in ps ], [ x[2] for x in ps ]
        for c1, c2 in zip(added_columns, valid):
            assert len(group) == len(c1), f"DataFrame length mismatch: {len(group)} != {len(c1)}\n\n{group}\n\n{c1}"
            assert len(group) == len(c2), f"ValidArray length mismatch: {len(group)} != {len(c2)}\n\n{group}\n\n{c2}"
        chunk = df_hstack([ group, *added_columns ])
        if args.filter:
            valid = functools.reduce(np.logical_and, valid)
            chunk = chunk.filter(valid)
        chunk = postfilter(chunk)
        if len(chunk) > 0:
            result_chunks.append(chunk)

    validate_missing_columns(result_chunks)
    df = pl.concat(result_chunks)
    df = postfilter(df)
    df = df.sort('pdbid', 'model', 'chain1', 'nr1', 'chain2', 'nr2')
    return df

def main(pool: Union[Pool, MockPool], args):
    if args.reference_basepairs:
        print("Loading metadata from ", args.reference_basepairs)
        ideal_basepairs = load_ideal_pairs(pool, args.reference_basepairs)
    else:
        ideal_basepairs = None

    if args.partition_input == 0:
        df = main_partition(pool, args, '', ideal_basepairs)
        save_output(args, df.lazy())
    else:
        pdbids = load_input_pdbids(args)
        partitions = {}
        for pdbid in pdbids:
            x = pdbid[:args.partition_input]
            partitions[x] = partitions.get(x, 0) + 1
        print(f"Partitioning {len(pdbids)} PDBIDs into {len(partitions)} partitions: {partitions}")
        p_results = dict()
        try:
            for p in sorted(partitions.keys()):
                print("Processing partition", p)
                p_result = main_partition(pool, args, p, ideal_basepairs)
                if len(p_result) > 0:
                    file = f"{args.output}_part{p}"
                    p_result.write_parquet(file)
                    p_results[p] = file

            all_df = pl.concat([ pl.scan_parquet(p, cache=False, low_memory=True) for p in p_results.values() ])
            save_output(args, all_df)
        finally:
            for p in p_results.values():
                os.remove(p)


def save_output(args, df: pl.LazyFrame):
    if args.output.endswith(".parquet"):
        df.sink_parquet(args.output)
    elif args.output != "/dev/null":
        df.sink_csv(args.output if args.output.endswith(".csv") else args.output + ".csv")
        df.sink_parquet(args.output + ".parquet")
    return df

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="""
        Adds geometric information to the specified basepairing CSV/Parquet/FR3D files.
        Added columns:
            * hb_0_length, hb_0_acceptor_angle, hb_0_donor_angle, ... - hydrogen bond lengths and angles between heavy atoms on both sides
            * hb_0_OOPA1, hb_0_OOPA2, ... - hydrogen bond / base plane angles
            * coplanarity_angle - angle between normals of the two base planes
            * coplanarity_shift1, coplanarity_shift2 - (minimum) distance between the pairing edge and the other base plane
            * coplanarity_edge_angle1 - angle between the pairing edge and the other base plane
            * C1_C1_distance - distance between the C1' atoms
            * C1_C1_total_angle - angle between the C1'-N bonds
            * Relative positions of N - C1' reference frames
                - x1, y1, z1, x2, y2, z2 - x, y, z coordinates of the N-C1' reference frame (1: left to right, 2: right to left)
                - yaw1, pitch1, roll1, yaw2, pitch2, roll2 - yaw, pitch, roll angles of the N-C1' reference frame
                - euler_phi, euler_theta, euler_psi - euler angles between the N-C1' reference frames (left to right, the other direction is  symmetric)
            * Distance to ideal basepair (if --reference-basepairs is specified)
                - rmsd_edge1, rmsd_edge2 - RMSD of the pairing edge, when fitted on the other bases
                - rmsd_C1N_frames1, rmsd_C1N_frames2 - RMSD of the C1', N1/N9, C6/C8 atoms, when fitted on the other bases
                - rmsd_all_base - "classic" RMSD of all atoms
        When --dssr-binary is specified, DSSR --analyze is executed to gain additional information:
            * dssr_pairing_type - pairing type according to DSSR (e.g. WC, Platform, ~rHoogsteen)
            * pairing_type, shear, stretch, stagger, buckle, propeller, opening, shift, slide, rise, tilt, roll, twist
        """,
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("inputs", nargs="+")
    parser.add_argument("--pdbcache", nargs="+", help="Directories to search for PDB files in order to avoid downloading. Last directory will be written to, if the structure is not found and has to be downloaded from RCSB. Also can be specified as PDB_CACHE_DIR env variable.")
    parser.add_argument("--output", "-o", required=True, help="Output CSV/Parquet file name. Both CSV and Parquet are always written.")
    parser.add_argument("--threads", type=int, default=1, help="Number of worker processes to spawn. Does not affect Polars threading, at the start, the process might use more threads than specified.")
    parser.add_argument("--export-only", default=False, action="store_true", help="Only re-export the input as CSV+Parquet files, do not calculate anything")
    parser.add_argument("--partition_input", default=0, type=int, help="Partition the input by N first characters of the PDBID. Reduces memory usage for large datasets.")
    parser.add_argument("--metadata", type=bool, default=True, help="Add deposition_date, resolution and structure_method columns")
    parser.add_argument("--dssr-binary", type=str, help="If specified, DSSR --analyze will be invoked for each structure and its results stored as 'dssr_' prefixed columns")
    parser.add_argument("--filter", default=False, action="store_true", help="Filter out rows for which the values could not be calculated")
    parser.add_argument("--postfilter-hb", default=None, type=float, help="Only include rows with at least 1 hydrogen bond length < X Å")
    parser.add_argument("--postfilter-shift", default=None, type=float, help="Only include rows with abs(coplanarity_shift) < X Å (both left and right)")
    parser.add_argument("--dedupe", default=False, action="store_true", help="Remove duplicate pairs, keep the one with preferred family (W > H > S), FR3D ordering, preferred base order (A > G > C > U), shorter bonds, or lower chain1,nr1")
    parser.add_argument("--reference-basepairs", type=str, help="output of pair_distributions.py with the 'nicest' basepairs, will be used as reference for RMSD calculation (rmsd_edge1, rmsd_edge2, rmsd_C1N_frames1, rmsd_C1N_frames2, rmsd_edge_C1N_frame, rmsd_all_base columns)")
    parser.add_argument("--disable-cross-symmetry", default=False, action="store_true", help="Skip basepairs with bases in different Asymmetrical Units")
    parser.add_argument("--override-pair-family", type=str, required=False, help="Ignore the pair family from the input and assume the specified one instead.")
    args = parser.parse_args()

    for x in args.pdbcache or []:
        pdb_utils.pdb_cache_dirs.append(os.path.abspath(x))
    os.environ["PDB_CACHE_DIR"] = ';'.join(pdb_utils.pdb_cache_dirs)

    multiprocessing.set_start_method("spawn")
    if args.threads == 1:
        main(MockPool(), args)
    else:
        with multiprocessing.Pool(processes=args.threads) as pool:
            main(pool, args)
