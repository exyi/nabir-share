import functools
import itertools
from multiprocessing.pool import Pool
from typing import Any, Callable, List, Optional, Tuple, TypeVar, Union
import Bio.PDB, Bio.PDB.Structure, Bio.PDB.Model, Bio.PDB.Residue, Bio.PDB.Atom, Bio.PDB.Chain, Bio.PDB.Entity
import os, sys, io, gzip, math, numpy as np
import polars as pl
from dataclasses import dataclass
import dataclasses
from requests import HTTPError
import multiprocessing
from pair_csv_parse import scan_pair_csvs
import pdb_utils
import pair_defs as pair_defs
from async_utils import MockPool
pdef = pair_defs

_sentinel = object()


@dataclass
class AltResidue:
    res: Bio.PDB.Residue.Residue
    alt: str

    @property
    def resname(self):
        return self.res.resname

    def transform_atom(self, a: Bio.PDB.Atom.DisorderedAtom) -> Any:
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
        if name in self.res:
            return self.transform_atom(self.res[name])
        else:
            if default is _sentinel:
                raise KeyError(f"Atom {name} not found in residue {self.res.full_id}")
            else:
                return default

def get_base_atoms(res: AltResidue) -> Optional[Tuple[Bio.PDB.Atom.Atom, List[Bio.PDB.Atom.Atom]]]:
    c1 = res.get_atom("C1'", None)
    if c1 is not None:
        planar_base_atoms = [ a for a in res.get_atoms() if not a.name.endswith("'") and a.element != "H" and a.element != "D" ]
        assert len(planar_base_atoms) >= 6
        return c1, planar_base_atoms
    return None

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
    arr = np.array([ a.coord for a in atoms ]) - c1.coord
    c1_bonded = np.linalg.norm(arr, axis=1) <= 1.6
    if not np.any(c1_bonded):
        raise ResideTransformError(f"No bonding atom found at {c1.full_id}")
    adj_matrix = np.linalg.norm(arr[:, None, :] - arr[None, :, :], axis=2) <= 1.6
    connected_component = bfs(adj_matrix, c1_bonded)

    filter = connected_component
    return [ a for a, b in zip(atoms, filter) if b ]

def normalize(v):
    return v / np.linalg.norm(v)

def orthonormal_basis(B):
    Q, R = np.linalg.qr(B)
    return Q

def get_reasonable_nucleotides(model: Bio.PDB.Model.Model) -> List[Bio.PDB.Residue.Residue]:
    nucleotides = []
    for r in model.get_residues():
        r: Bio.PDB.Residue.Residue = r
        if "C1'" not in r:
            continue # TODO: is this naming consistent in all structures?
        if r.resname == "HOH":
            continue
        atoms = [ a for a in r.get_atoms() if a.element != "H" and a.element != "D" ]
        if len(atoms) < 6:
            continue
        if sum(1 for a in atoms if a.element == "C") < 4:
            continue
        if sum(1 for a in atoms if a.element == "N") < 2:
            continue
        
    return nucleotides

@dataclass
class TranslationThenRotation:
    translation: np.ndarray
    rotation: np.ndarray
    def __post_init__(self):
        assert self.rotation.shape == (3, 3)
        assert self.translation.shape == (3,)

@dataclass
class ResiduePosition:
    # translation of C1 to origin
    origin: np.ndarray
    # rotation around C1 to align C1'-N bond to x-axis and the plane to X/Y axes
    rotation: np.ndarray
    plane_basis: np.ndarray
    plane_normal_vector: np.ndarray
    projection_matrix: np.ndarray
    c1_pos: Optional[np.ndarray] = None
    n_pos: Optional[np.ndarray] = None

    def __post_init__(self):
        assert self.rotation.shape == (3, 3)
        assert self.plane_basis.shape == (3, 2)

    def plane_projection(self, point: np.ndarray) -> np.ndarray:
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
        normal_angle = math.degrees(math.acos(np.dot(normalize(line), self.plane_normal_vector)))
        return 90 - normal_angle

def fit_plane_magic_svd(atoms: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """
    Returns centroid + S (= basis1, basis2, normal)
    https://math.stackexchange.com/a/99317
    """
    assert atoms.shape[1] == 3
    svd = np.linalg.svd((atoms - np.mean(atoms, axis=0, keepdims=True)).T)
    left_singular_vectors = svd[0]
    assert left_singular_vectors.shape == (3, 3)
    return np.mean(atoms, axis=0), left_singular_vectors

def get_residue_posinfo(res: AltResidue) -> ResiduePosition:
    """
    Returns a tuple of (translation, rotation) that moves the C1' atom to the origin and the C1'-N bond to the x-axis and all planar atoms with (near) 0 z-coordinate.
    """
    x = get_base_atoms(res)
    if x is None:
        raise ResideTransformError(f"No C1' atom found at {res.res.full_id}")
    c1, planar_atoms = x
    nX = res.get_atom("N1", None) or res.get_atom("N9", None)
    planar_atoms = filter_atoms(c1, planar_atoms)
    
    atom_names = np.array([ a.name for a in planar_atoms ])
    atom_elements = np.array([ a.element for a in planar_atoms ])
    atoms = np.array([ a.coord for a in planar_atoms ])

    # dist_matrix = np.linalg.norm(atoms[:, None, :] - atoms[None, :, :], axis=2)
    # dist_matrix[np.arange(len(atoms)), np.arange(len(atoms))] = 1000_000

    # check enough distance from C1'
    # if np.any(np.linalg.norm(atoms, axis=1) < 1.3):
    #     raise ResideTransformError(f"Atoms too close to origin at {res.res.full_id}")
    
    # * fit a plane through the atoms
    # fitted_plane, sum_error, _, _ = np.linalg.lstsq(np.concatenate([atoms[:, 1:], np.ones((len(planar_atoms), 1))], axis=1), -atoms[:, 0], rcond=None) # x + ay + bz + c = 0
    # fitted_plane, sum_error, _, _ = np.linalg.lstsq(np.concatenate([atoms, np.ones((atoms.shape[0], 1))], axis=1), np.zeros(atoms.shape[0]), rcond=None) # ax + by + cz + d = 0
    # fitted_plane, sum_error, _, _ = np.linalg.lstsq(atoms[:, :2], -atoms[:, 2], rcond=None)
    # fitted_plane, sum_error, _, _ = np.linalg.lstsq(np.concatenate([atoms, np.ones, np.zeros(len(planar_atoms)), rcond=None) # ax + by + cz = 0

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
        # rot_matrix1 = np.concatenate([plane_basis, normal_v.reshape(3, 1)], axis=1)
        assert np.sum((atoms_c @ rot_matrix1)[:, 2] ** 2) < 1.3 * sum_error
        
        assert np.allclose(rot_matrix1 @ rot_matrix1.T, np.eye(3), atol=1e-5)
        assert np.allclose(np.linalg.det(rot_matrix1), 1, atol=1e-5)

        atoms_c_rot1 = atoms_c @ rot_matrix1

        assert np.allclose(np.linalg.norm(atoms_c, axis=1), np.linalg.norm(atoms_c_rot1, axis=1), atol=1e-5)

        # orient C1'-N bond to x-axis (using 2D rotation on X,Y)
        c1_c_rot = (c1.coord - origin) @ rot_matrix1
        x_vector = -normalize(c1_c_rot)
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
    Orients residue so that the C1-N bond is aligned with the x-axis, the N1-C2 or N9-C8 bond is towards Y axis and remaining is Z
    """
    c1 = res.get_atom("C1'", None)
    if res.get_atom("N9", None) is not None:
        n = res.get_atom("N9", None)
        c2 = res.get_atom("C8", None)
    else:
        n = res.get_atom("N1", None)
        c2 = res.get_atom("C6", None) # ale tady mi copilot dal hajzl C dvojku
    if c1 is None or n is None or c2 is None:
        raise ResideTransformError(f"Missing atoms in residue {res.res.full_id}")
    translation = -n.coord # tvl copilot toto dal asi, cool priklad do appendix AI
    x = c1.coord - n.coord
    x /= np.linalg.norm(x)
    y = c2.coord - n.coord
    y -= np.dot(y, x) * x
    y /= np.linalg.norm(y)
    z = np.cross(y, x)
    z /= np.linalg.norm(z)
    rotation = np.array([x, y, z]).T

    test = AltResidue(transform_residue(res.res, translation, rotation), res.alt)
    assert np.all(test.get_atom(n.name).coord < 0.001)
    assert np.all(test.get_atom(c1.name).coord[1:] < 0.001)
    assert 1.2 < test.get_atom(c1.name).coord[0] < 1.9, f"Unexpected C1' coordinates (x is bad): {test.get_atom(c1.name).coord}"
    assert np.all(test.get_atom(c2.name).coord[2] < 0.001)

    return TranslationThenRotation(translation, rotation)

def get_C1_N_yaw_pitch_roll(rot1: np.ndarray, rot2: np.ndarray) -> Tuple[float, float, float]:
    matrix = rot1.T @ rot2
    # https://stackoverflow.com/a/37558238
    yaw = math.degrees(math.atan2(matrix[1, 0], matrix[0, 0]))
    pitch = math.degrees(math.atan2(-matrix[2, 0], math.sqrt(matrix[2, 1]**2 + matrix[2, 2]**2)))
    roll = math.degrees(math.atan2(matrix[2, 1], matrix[2, 2]))
    return yaw, pitch, roll
def get_C1_N_euler_angles(rot1: np.ndarray, rot2: np.ndarray) -> Tuple[float, float, float]:
    matrix = rot1.T @ rot2
    # https://en.wikipedia.org/wiki/Rotation_formalisms_in_three_dimensions#Rotation_matrix_%E2%86%92_Euler_angles_(z-x-z_extrinsic)
    eu1 = math.degrees(math.atan2(matrix[2, 0], matrix[2, 1]))
    eu2 = math.degrees(math.acos(matrix[2, 2]))
    eu3 = -math.degrees(math.atan2(matrix[0, 2], matrix[1, 2]))
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

def transform_residue(r: Bio.PDB.Residue.Residue, translation: np.ndarray, rotation: np.ndarray):
    r = r.copy()
    r.transform(np.eye(3), translation)
    r.transform(rotation, np.zeros(3))
    return r

def to_float32_df(df: pl.DataFrame) -> pl.DataFrame:
    return df.with_columns([
        pl.col(c).cast(pl.Float32).alias(c) for c in df.columns if df[c].dtype == pl.Float64
    ])

@dataclass
class PairStats:
    # buckle: float
    coplanarity_angle: Optional[float] = None
    coplanarity_shift1: Optional[float] = None
    coplanarity_shift2: Optional[float] = None
    coplanarity_edge_angle1: Optional[float] = None
    coplanarity_edge_angle2: Optional[float] = None
    C1_C1_distance: Optional[float] = None
    C1_C1_total_angle: Optional[float] = None
    C1_C1_yaw1: Optional[float] = None
    C1_C1_pitch1: Optional[float] = None
    C1_C1_roll1: Optional[float] = None
    C1_C1_yaw2: Optional[float] = None
    C1_C1_pitch2: Optional[float] = None
    C1_C1_roll2: Optional[float] = None
    C1_C1_euler_phi: Optional[float] = None
    C1_C1_euler_theta: Optional[float] = None
    C1_C1_euler_psi: Optional[float] = None
    # opening: float
    # nearest_distance: float

def _linalg_norm(v) -> float:
    return float(np.linalg.norm(v)) # typechecker se asi posral nebo fakt nevím z jakých halucinací usoudil, že výsledek np.linalg.norm nejde násobit, porovnávat nebo nic

def _linalg_normalize(v) -> np.ndarray:
    return v / np.linalg.norm(v)

def maybe_agg(agg, array):
    if len(array) == 0:
        return None
    return agg(array)

def abs_min(array):
    ix = np.argmin(np.abs(array))
    return array[ix]

def get_edge_atoms(res: AltResidue, edge_name: str) -> List[Bio.PDB.Atom.Atom]:
    edge = pair_defs.base_edges.get(resname_map.get(res.resname, res.resname), {}).get(edge_name, [])
    return [ a for aname in edge if (a:=res.get_atom(aname, None)) is not None ]

def get_edge_edges(edge: List[Bio.PDB.Atom.Atom]) -> Optional[Tuple[Bio.PDB.Atom.Atom, Bio.PDB.Atom.Atom]]:
    if not edge or len(edge) < 2:
        return None
    polar_atoms = [ a for a in edge if a.element not in ("H", "C") ]
    if len(polar_atoms) < 2:
        return (edge[0], edge[-1])
    return (polar_atoms[0], polar_atoms[-1])

def calc_pair_stats(pair_type: pair_defs.PairType, res1: AltResidue, res2: AltResidue) -> tuple[Optional[PairStats], Optional[ResiduePosition], Optional[ResiduePosition]]:
    p1 = try_get_residue_posinfo(res1)
    p2 = try_get_residue_posinfo(res2)
    rot_trans1 = try_x(get_residue_posinfo_C1_N, res1)
    rot_trans2 = try_x(get_residue_posinfo_C1_N, res2)

    if p1 is None or p2 is None:
        return None, p1, p2
    res1_ = AltResidue(transform_residue(res1.res, -p1.origin, p1.rotation), res1.alt)
    res2_ = AltResidue(transform_residue(res2.res, -p1.origin, p1.rotation), res2.alt)
    edge1 = get_edge_atoms(res1, pair_type.edge1_name)
    edge1_edges = get_edge_edges(edge1)
    edge2 = get_edge_atoms(res2, pair_type.edge2_name)
    edge2_edges = get_edge_edges(edge2)

    out = PairStats()

    out.coplanarity_angle = math.degrees(np.arccos(np.dot(p1.plane_normal_vector, p2.plane_normal_vector)))
    
    out.coplanarity_shift1 = maybe_agg(abs_min, [p2.plane_get_deviation(a.coord) for a in edge1])
    out.coplanarity_shift2 = maybe_agg(abs_min, [p1.plane_get_deviation(a.coord) for a in edge2])
    out.coplanarity_edge_angle1 = p2.get_relative_line_rotation(*(a.coord for a in edge1_edges)) if edge1_edges is not None else None
    out.coplanarity_edge_angle2 = p1.get_relative_line_rotation(*(a.coord for a in edge2_edges)) if edge2_edges is not None else None

    out.C1_C1_distance = _linalg_norm(res1.get_atom("C1'").coord - res2.get_atom("C1'").coord)
    if p1.n_pos is not None and p2.n_pos is not None and p1.c1_pos is not None and p2.c1_pos is not None:
        out.C1_C1_total_angle = math.degrees(np.arccos(np.dot(_linalg_normalize(p1.n_pos - p1.c1_pos), _linalg_normalize(p2.c1_pos - p2.n_pos))))
    
    if rot_trans1 and rot_trans2:
        out.C1_C1_yaw1, out.C1_C1_pitch1, out.C1_C1_roll1 = get_C1_N_yaw_pitch_roll(rot_trans1.rotation, -rot_trans2.rotation)
        out.C1_C1_yaw2, out.C1_C1_pitch2, out.C1_C1_roll2 = get_C1_N_yaw_pitch_roll(rot_trans2.rotation, -rot_trans1.rotation)
        out.C1_C1_euler_phi, out.C1_C1_euler_theta, out.C1_C1_euler_psi = get_C1_N_euler_angles(rot_trans1.rotation, -rot_trans2.rotation)

    return out, p1, p2

def get_angle(atom1: Bio.PDB.Atom.Atom, atom2: Bio.PDB.Atom.Atom, atom3: Bio.PDB.Atom.Atom) -> float:
    v1 = atom1.coord - atom2.coord
    v2 = atom3.coord - atom2.coord
    return math.degrees(np.arccos(np.dot(v1, v2) / (_linalg_norm(v1) * _linalg_norm(v2))))

def get_distance(atom1: Bio.PDB.Atom.Atom, atom2: Bio.PDB.Atom.Atom) -> float:
    return _linalg_norm(atom1.coord - atom2.coord)


@dataclass
class HBondStats:
    length: float
    donor_angle: float
    acceptor_angle: float
    donor_OOPA: Optional[float] = None
    acceptor_OOPA: Optional[float] = None

def hbond_stats(
    atom0: Bio.PDB.Atom.Atom, atom1: Bio.PDB.Atom.Atom, atom2: Bio.PDB.Atom.Atom, atom3: Bio.PDB.Atom.Atom,
    donor_plane: Optional[ResiduePosition],
    acceptor_plane: Optional[ResiduePosition]
) -> HBondStats:
    assert _linalg_norm(atom0.coord - atom1.coord) <= 1.6, f"atoms 0,1 not bonded: {atom0.full_id} {atom1.full_id} ({np.linalg.norm(atom0.coord - atom1.coord)} > 1.6, {atom0.coord} {atom1.coord})"
    assert _linalg_norm(atom2.coord - atom3.coord) <= 1.6, f"atoms 2,3 not bonded: {atom2.full_id} {atom3.full_id} ({np.linalg.norm(atom0.coord - atom1.coord)} > 1.6, {atom0.coord} {atom1.coord})"
    assert _linalg_norm(atom1.coord - atom2.coord) > 1.55, f"atoms too close for h-bond: {atom1.full_id} {atom2.full_id} ({np.linalg.norm(atom1.coord - atom2.coord)} < 2, {atom1.coord} {atom2.coord})"
    assert _linalg_norm(atom1.coord - atom2.coord) < 10, f"atoms too far for h-bond: ({np.linalg.norm(atom1.coord - atom2.coord)} > 6, {atom1.coord} {atom2.coord})"


    result = HBondStats(
        length=get_distance(atom1, atom2),
        donor_angle=get_angle(atom0, atom1, atom2),
        acceptor_angle=get_angle(atom3, atom2, atom1)
    )

    if donor_plane is not None:
        result.donor_OOPA = donor_plane.get_relative_line_rotation(atom1.coord, atom2.coord)
    if acceptor_plane is not None:
        result.acceptor_OOPA = acceptor_plane.get_relative_line_rotation(atom2.coord, atom1.coord)

    return result

resname_map = {
    'DT': 'U',
    'DC': 'C',
    'DA': 'A',
    'DG': 'G',
    'DU': 'U',
    'T': 'U',
}

def get_hbond_stats(pair_type: pair_defs.PairType, r1: AltResidue, r2: AltResidue, rp1: Optional[ResiduePosition], rp2: Optional[ResiduePosition]) -> Optional[List[Optional[HBondStats]]]:
    # if r1.resname > r2.resname:
    #     r1, r2 = r2, r1
    hbonds = pdef.get_hbonds(pair_type, throw=False)
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
        hbond_stats(atoms[0], atoms[1], atoms[2], atoms[3], rp1 if hb[1][0] == 'A' else rp2, rp1 if hb[2][0] == 'A' else rp2) if None not in atoms else None
        for atoms, hb in zip(bonds_atoms, hbonds)
    ]

def get_atom_df(model: Bio.PDB.Model.Model):
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

def get_residue(structure: Bio.PDB.Structure.Structure, model, chain, res, ins, alt) -> AltResidue:
    res = int(res)
    ins = ins or ' '
    ch: Bio.PDB.Chain.Chain = structure[model-1][chain]
    found = ch.child_dict.get((' ', res, ins), None)
    if found is None:
        found = next((r for r in ch.get_residues() if r.id[1] == res and r.id[2] == ins), None)
    if found is None:
        raise KeyError(f"{chain}.{res}{ins}")
    return AltResidue(found, alt)


def get_stats_for_csv(df_: pl.DataFrame, structure: Bio.PDB.Structure.Structure, global_pair_type: Optional[str]):
    df = df_.with_row_count()
    pair_type_col = df["type"] if global_pair_type is None else itertools.repeat(global_pair_type)
    for (pair_family, i, pdbid, model, chain1, res1, ins1, alt1, chain2, res2, ins2, alt2) in zip(pair_type_col, df["row_nr"], df["pdbid"], df["model"], df["chain1"], df["nr1"], df["ins1"], df["alt1"], df["chain2"], df["nr2"], df["ins2"], df["alt2"]):
        if (chain1, res1, ins1, alt1) == (chain2, res2, ins2, alt2):
            continue
        try:
            assert structure.id.lower() == pdbid.lower(), f"pdbid mismatch: {structure.id} != {pdbid}"
            ins1 = ins1.strip()
            ins2 = ins2.strip()

            try:
                r1 = get_residue(structure, model, chain1, res1, ins1, alt1)
            except KeyError as keyerror:
                print(f"Could not find residue1 {pdbid}.{model}.{chain1}.{str(res1)+ins1} ({keyerror=})")
                continue
            try:
                r2 = get_residue(structure, model, chain2, res2, ins2, alt2)
            except KeyError as keyerror:
                print(f"Could not find residue2 {pdbid}.{model}.{chain2}.{str(res2)+ins2} in {pdbid} ({keyerror=})")
                continue
            pair_type = pair_defs.PairType.create(pair_family, r1.resname, r2.resname, name_map=resname_map)
            stats, rp1, rp2 = None, None, None
            stats, rp1, rp2 = calc_pair_stats(pair_type, r1, r2)
            hbonds = get_hbond_stats(pair_type, r1, r2, rp1, rp2)
            if hbonds is not None:
                yield i, hbonds, stats
        except AssertionError as e:
            print(f"Assertion error {e} on row:\n{to_csv_row(df_, i, 15)}")
            continue
        except Exception as keyerror:
            print(f"Error on row:\n{to_csv_row(df_, i, 15)}")
            import traceback
            print(keyerror)
            print(traceback.format_exc())
            print("Continuing...")
            continue

def remove_duplicate_pairs(df: pl.DataFrame):
    def pair_id(type, res1, res2, chain1, nr1, ins1, alt1, chain2, nr2, ins2, alt2):
        pair = [ (chain1, nr1, ins1, alt1), (chain2, nr2, ins2, alt2) ]
        pair.sort()
        bases = resname_map.get(res1, res1) + "-" + resname_map.get(res2, res2)
        return type + "-" + bases + "|".join(tuple(str(x) for x in (pair[0] + pair[1])))
    pair_ids = [ pair_id(type, res1, res2, chain1, nr1, ins1, alt1, chain2, nr2, ins2, alt2) for type, res1, res2, chain1, nr1, ins1, alt1, chain2, nr2, ins2, alt2 in zip(df["type"], df["res1"], df["res2"], df["chain1"], df["nr1"], df["ins1"], df["alt1"], df["chain2"], df["nr2"], df["ins2"], df["alt2"]) ]
    df = df.with_columns(
        pl.Series(pair_ids, dtype=pl.Utf8).alias("_tmp_pair_id")
    ).with_row_count("_tmp_row_nr")
    score = pl.lit(0, dtype=pl.Float64)
    for col in df.columns:
        if col.startswith("hb_") and col.endswith("_length"):
            score += pl.col(col).fill_null(100)
        if col.startswith("dssr_"):
            score += pl.col(col).is_null().cast(pl.Float64) * 0.03
        if col == "coplanarity_angle":
            score += pl.col(col).is_null().cast(pl.Float64)

    df = df.sort([score, "chain1", "nr1", "ins1", "alt1", "chain2", "nr2", "ins2", "alt2"])
    df = df.unique(["type", "pdbid", "model", "_tmp_pair_id"], keep="first", maintain_order=True)
    df = df.sort("_tmp_row_nr")
    df = df.drop(["_tmp_pair_id", "_tmp_row_nr"])
    return df

def get_max_bond_count(df: pl.DataFrame):
    pair_types = list(set(pair_defs.PairType.create(type, res1, res2, name_map=resname_map) for type, res1, res2 in df[["type", "res1", "res2"]].unique().iter_rows()))
    bond_count = max(len(pair_defs.get_hbonds(p, throw=False)) for p in pair_types)
    print("Analyzing pair types:", *pair_types, "with bond count =", bond_count)
    return max(3, bond_count)

def backbone_columns(df: pl.DataFrame, structure: Optional[Bio.PDB.Structure.Structure]) -> list[pl.Series]:
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
        if (nr1, ins1 or ' ') in chain_index[(model, chain1)] and (nr2, ins2 or ' ') in chain_index[(model, chain2)]
    )

    for i, (pdbid, model, chain1, nr1, ins1, alt1, chain2, nr2, ins2, alt2) in enumerate(zip(df["pdbid"], df["model"], df["chain1"], df["nr1"], df["ins1"], df["alt1"], df["chain2"], df["nr2"], df["ins2"], df["alt2"])):

        res_ix1 = chain_index[(model, chain1)].get((nr1, ins1 or ' '), None)
        res_ix2 = chain_index[(model, chain2)].get((nr2, ins2 or ' '), None)
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

        if antipar and par:
            print("Interesting pair:", pdbid, model, chain1, nr1, ins1, chain2, nr2, ins2, "both parallel and antiparallel")
    return [ is_dinucleotide, is_parallel ]



def export_stats_csv(pdbid, df: pl.DataFrame, add_metadata_columns: bool, pair_type: Optional[str], max_bond_count: int) -> Tuple[str, pl.DataFrame, np.ndarray]:
    bond_params = [ x.name for x in dataclasses.fields(HBondStats) ]
    valid = np.zeros(len(df), dtype=np.bool_)
    columns: list[pl.Series] = [
        pl.Series(f"hb_{i}_{p}", [ None ] * len(df), dtype=pl.Float32)
        for i in range(max_bond_count)
        for p in bond_params ]
    pair_stats: list[Optional[PairStats]] = [ None ] * len(df)

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

    if structure is not None:
        for i, hbonds, stats in get_stats_for_csv(df, structure, pair_type):
            valid[i] = True
            for j, s in enumerate(hbonds):
                if s is None:
                    continue

                for param_i, param in enumerate(bond_params):
                    columns[j * len(bond_params) + param_i][i] = getattr(s, param)

            pair_stats[i] = stats

    result_df = pl.DataFrame(columns).hstack(pl.DataFrame([
            p or PairStats() for p in pair_stats
        ]))
    result_df = to_float32_df(result_df)
    if add_metadata_columns:
        h = structure.header if structure is not None else dict()
        result_df = result_df.with_columns(
            pl.lit(h.get('deposition_date', None), dtype=pl.Utf8).alias("deposition_date"),
            pl.lit(h.get('name', None), dtype=pl.Utf8).alias("structure_name"),
            pl.lit(h.get('structure_method', None), dtype=pl.Utf8).alias("structure_method"),
            pl.lit(h.get('resolution', None), dtype=pl.Float32).alias("resolution")
        )
    
    result_df = result_df.with_columns(backbone_columns(df, structure))
    assert len(result_df) == len(df)
    return pdbid, result_df, valid


def df_hstack(columns: list[pl.DataFrame]) -> pl.DataFrame:
    return functools.reduce(pl.DataFrame.hstack, columns)

def load_inputs(pool: Union[Pool, MockPool], args) -> pl.DataFrame:
    inputs: list[str] = args.inputs
    if len(inputs) == 0:
        raise ValueError("No input files specified")
    
    if inputs[0].endswith(".parquet") or inputs[0].endswith('.csv'):
        print(f'Loading basepairing CSV files')
        df = scan_pair_csvs(args.inputs).sort('pdbid', 'model', 'nr1', 'nr2').collect()
    elif inputs[0].endswith("_basepair.txt") or inputs[0].endswith("_basepair_detail.txt"):
        print(f"Loading {len(inputs)} basepair files")
        import fr3d_parser
        df = fr3d_parser.read_fr3d_files_df(pool, args.inputs, filter=pl.col("symmetry_operation") == '').sort('pdbid', 'model', 'nr1', 'nr2')
    else:
        raise ValueError("Unknown input file type")
    if "type" not in df.columns:
        if not args.pair_type:
            raise ValueError("Input does not contain type column and --pair-type was not specified")
        df = df.select(pl.lit(args.pair_type).alias("type"), pl.col("*"))
    return df

def validate_missing_columns(chunks: list[pl.DataFrame]):
    all_columns = set(col for chunk in chunks for col in chunk.columns)
    for c in chunks:
        if not all_columns.issubset(c.columns):
            print(f"Chunk with {set(c['pdbid'].to_numpy())} is missing columns: {all_columns.difference(set(c.columns))}")

def main(pool: Union[Pool, MockPool], args):
    df = load_inputs(pool, args)
    max_bond_count = get_max_bond_count(df)
    groups = list(df.group_by(pl.col("pdbid")))

    processes = []

    processes.append([
        pool.apply_async(export_stats_csv, args=[pdbid, group, args.metadata, args.pair_type, max_bond_count])
        for (pdbid,), group in groups # type: ignore
        # for chunk in group.iter_slices(n_rows=100)
    ])
    if args.dssr_binary is not None:
        import dssr_wrapper
        processes.append([
            pool.apply_async(dssr_wrapper.add_dssr_info, args=[pdbid, group, args.dssr_binary])
            for ((pdbid,), group) in groups # type: ignore
        ])

    result_chunks = []
    for ((_pdbid,), group), p in zip(groups, zip(*processes)):
        p = [ x.get() for x in p ]
        pdbids = [ x[0] for x in p ]
        assert pdbids == [ _pdbid ] * len(pdbids), f"pdbid mismatch: {_pdbid} x {pdbids}"
        added_columns, valid = [ x[1] for x in p ], [ x[2] for x in p ]
        for c1, c2 in zip(added_columns, valid):
            assert len(group) == len(c1), f"DataFrame length mismatch: {len(group)} != {len(c1)}\n\n{group}\n\n{c1}"
            assert len(group) == len(c2), f"ValidArray length mismatch: {len(group)} != {len(c2)}\n\n{group}\n\n{c2}"
        chunk = df_hstack([ group, *added_columns ])
        if args.filter:
            valid = functools.reduce(np.logical_and, valid)
            chunk = chunk.filter(valid)
        if len(chunk) > 0:
            result_chunks.append(chunk)

    validate_missing_columns(result_chunks)
    df = pl.concat(result_chunks)
    if args.dedupe:
        df = remove_duplicate_pairs(df)
    df = df.sort('pdbid', 'model', 'nr1', 'nr2')
    if args != "/dev/null":
        df.write_csv(args.output if args.output.endswith(".csv") else args.output + ".csv")
        df.write_parquet(args.output + ".parquet")
    return df

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="""
        Adds geometric information to the specified pairing CSV files.
        Added columns:
            * hb_0_length, hb_0_acceptor_angle, hb_0_donor_angle, ... - hydrogen bond lengths and angles between heavy atoms on both sides
            * coplanarity_angle - angle between planes of the two nucleotides
        When --dssr-binary is specified, DSSR --analyze is executed to gain additional information:
            * dssr_pairing_type - pairing type according to DSSR (e.g. WC, Platform, ~rHoogsteen)
            * pairing_type, shear, stretch, stagger, buckle, propeller, opening, shift, slide, rise, tilt, roll, twist
        """,
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("inputs", nargs="+")
    parser.add_argument("--pdbcache", nargs="+", help="Directories to search for PDB files in order to avoid downloading. Last directory will be written to, if the structure is not found and has to be downloaded.")
    parser.add_argument("--output", "-o", required=True, help="Output CSV/Parquet file name")
    parser.add_argument("--threads", type=int, default=1, help="Maximum parallelism - number of worker processes to spawn")
    parser.add_argument("--metadata", type=bool, default=True, help="Add deposition_date, resolution and structure_method columns")
    parser.add_argument("--dssr-binary", type=str, help="If specified, DSSR --analyze will be invoked for each structure and its results stored as 'dssr_' prefixed columns")
    parser.add_argument("--filter", default=False, action="store_true", help="Filter out rows for which the values could not be calculated")
    parser.add_argument("--dedupe", default=False, action="store_true", help="Remove duplicate pairs, keep the one with shorter bonds or lower chain1,nr1")
    parser.add_argument("--pair-type", type=str, required=False)
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
