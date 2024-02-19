import functools
import itertools
from multiprocessing.pool import Pool
from typing import Any, List, Optional, Tuple, Union
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
            # if self.alt == '':
            #     return a.disordered_get('A')
            if not a.disordered_has_id(self.alt):
                raise KeyError(f"Atom {self.res.full_id}:{a.id} is disordered, but alt='{self.alt}' must be one of {a.disordered_get_id_list()}")
            return a.disordered_get(self.alt)
        else:
            return a
    def get_atoms(self):
        for a in self.res.get_atoms():
            yield self.transform_atom(a)
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
class ResiduePosition:
    # translation of C1 to origin
    translation: np.ndarray
    # rotation around C1 to align C1'-N bond to x-axis and the plane to X/Y axes
    rotation: np.ndarray
    plane_basis: np.ndarray
    plane_normal_vector: np.ndarray

    def __post_init__(self):
        assert self.rotation.shape == (3, 3)
        assert self.plane_basis.shape == (3, 2)

    def plane_projection(self, point: np.ndarray) -> np.ndarray:
        assert point.shape[-1] == 3
        tpoint = point + self.translation
        return np.dot(tpoint, self.plane_basis[:, 0]) * self.plane_basis[:, 0] + np.dot(tpoint, self.plane_basis[:, 1]) * self.plane_basis[:, 1] - self.translation

    def plane_get_deviation(self, point: np.ndarray) -> float:
        """
        Returns the absolute distance of the point from the plane (in Å).
        """
        assert point.shape == (3,)
        return _linalg_norm(point - self.plane_projection(point))
    
    def get_projection_squish_angle(self, point1: np.ndarray, point2: np.ndarray) -> float:
        """
        Returns the ratio of the distance between the projections and the distance between the points.
        """
        assert point1.shape == (3,) == point2.shape
        proj1 = self.plane_projection(point1)
        proj2 = self.plane_projection(point2)
        dist_cmp = _linalg_norm(proj1 - proj2) / _linalg_norm(point1 - point2)
        return min(1, max(0, dist_cmp))


def get_residue_posinfo(res: AltResidue) -> ResiduePosition:
    """
    Returns a tuple of (translation, rotation) that moves the C1' atom to the origin and the C1'-N bond to the x-axis and all planar atoms with (near) 0 z-coordinate.
    """
    x = get_base_atoms(res)
    if x is None:
        raise ResideTransformError(f"No C1' atom found at {res.res.full_id}")
    c1, planar_atoms = x
    planar_atoms = filter_atoms(c1, planar_atoms)
    
    translation: np.ndarray = -c1.coord
    atom_names = np.array([ a.name for a in planar_atoms ])
    atom_elements = np.array([ a.element for a in planar_atoms ])
    atoms = np.array([ a.coord for a in planar_atoms ]) + translation

    # dist_matrix = np.linalg.norm(atoms[:, None, :] - atoms[None, :, :], axis=2)
    # dist_matrix[np.arange(len(atoms)), np.arange(len(atoms))] = 1000_000

    # find the bonding N atom (i.e. distance <= 1.6Å)
    if np.any(np.linalg.norm(atoms, axis=1) < 1.3):
        raise ResideTransformError(f"Atoms too close to origin at {res.res.full_id}")
    
    # * fit a plane through the atoms
    # fitted_plane, sum_error, _, _ = np.linalg.lstsq(np.concatenate([atoms[:, :2], np.ones((len(planar_atoms), 1))], axis=1), atoms[:, 2], rcond=None)
    fitted_plane, sum_error, _, _ = np.linalg.lstsq(atoms[:, :2], -atoms[:, 2], rcond=None)
    if sum_error > 0.5 * len(planar_atoms): # TODO: reasonable threshold?
        raise ResideTransformError(f"Residue really doesn't seem planar: {res.res.full_id}, plane RMSE = {np.sqrt(sum_error / len(planar_atoms))}, planar atoms = {list(atom_names)}")

    # rotate to align the fitted_plane to x, y
    plane_basis = orthonormal_basis(np.array([
        [1, 0, -fitted_plane[0]],
        [0, 1, -fitted_plane[1]]
    ]).T)
    # rotation matrix is an orthonormal matrix with det=1. Cross product gives us the remaining vector orthogonal to the plane.
    rot_matrix1 = np.concatenate([plane_basis, np.cross(plane_basis[:, 0], plane_basis[:, 1]).reshape(3, 1)], axis=1)
    
    assert np.allclose(rot_matrix1 @ rot_matrix1.T, np.eye(3))
    assert np.allclose(np.linalg.det(rot_matrix1), 1)

    atoms2 = atoms @ rot_matrix1

    assert np.allclose(np.linalg.norm(atoms, axis=1), np.linalg.norm(atoms2, axis=1))

    # orient C1'-N bond to x-axis (using 2D rotation on X,Y)
    c1_bonded = np.linalg.norm(atoms, axis=1) <= 1.6
    c1_bonded_atom = (list(atoms2[c1_bonded & (atom_elements == "N")]) or list(atoms2[c1_bonded]))[0]
    x_vector = normalize(c1_bonded_atom[0:2])
    rot_matrix2 = np.array([
        [x_vector[0], -x_vector[1], 0],
        [x_vector[1], x_vector[0],  0],
        [ 0,          0,            1]
    ])
    atoms3 = atoms2 @ rot_matrix2

    return ResiduePosition(
        translation,
        rot_matrix1 @ rot_matrix2,
        plane_basis,
        np.cross(plane_basis[:, 0], plane_basis[:, 1])
    )

def try_get_residue_posinfo(res: AltResidue, warning = None) -> Optional[ResiduePosition]:
    try:
        return get_residue_posinfo(res)
    except ResideTransformError as e:
        if warning is not None:
            warning(str(e))
        else:
            print(f"WARNING(residue = {res.res.full_id} {res.alt}): {e}")
        return None

def transform_residue(r: Bio.PDB.Residue.Residue, translation: np.ndarray, rotation: np.ndarray):
    r = r.copy()
    r.transform(np.eye(3), translation)
    r.transform(rotation, np.zeros(3))
    return r

@dataclass
class PairStats:
    # buckle: float
    coplanarity_angle: Optional[float]
    coplanarity_shift1: Optional[float]
    coplanarity_shift2: Optional[float]
    coplanarity_edge_angle1: Optional[float]
    coplanarity_edge_angle2: Optional[float]
    # opening: float
    # nearest_distance: float

def _linalg_norm(v) -> float:
    return float(np.linalg.norm(v)) # typechecker se asi posral nebo fakt nevím z jakých halucinací usoudil, že výsledek np.linalg.norm nejde násobit, porovnávat nebo nic

def maybe_min(array):
    if len(array) == 0:
        return None
    return min(array)

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

def calc_pair_stats(pair_type: pair_defs.PairType, res1: AltResidue, res2: AltResidue) -> Optional[PairStats]:
    p1 = try_get_residue_posinfo(res1)
    p2 = try_get_residue_posinfo(res2)
    if p1 is None or p2 is None:
        return None
    res1_ = AltResidue(transform_residue(res1.res, p1.translation, p1.rotation), res1.alt)
    res2_ = AltResidue(transform_residue(res2.res, p1.translation, p1.rotation), res2.alt)
    edge1 = get_edge_atoms(res1, pair_type.edge1_name)
    edge1_edges = get_edge_edges(edge1)
    edge2 = get_edge_atoms(res2, pair_type.edge2_name)
    edge2_edges = get_edge_edges(edge2)

    copl_angle = math.degrees(np.arccos(np.dot(p1.plane_normal_vector, p2.plane_normal_vector)))
    
    copl_shift1 = maybe_min([p2.plane_get_deviation(a.coord) for a in edge1])
    copl_shift2 = maybe_min([p1.plane_get_deviation(a.coord) for a in edge2])
    copl_edge_a1 = math.degrees(math.acos(p2.get_projection_squish_angle(*(a.coord for a in edge1_edges)))) if edge1_edges is not None else None
    copl_edge_a2 = math.degrees(math.acos(p1.get_projection_squish_angle(*(a.coord for a in edge2_edges)))) if edge2_edges is not None else None
    return PairStats(copl_angle, copl_shift1, copl_shift2, copl_edge_a1, copl_edge_a2)

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

def hbond_stats(atom0: Bio.PDB.Atom.Atom, atom1: Bio.PDB.Atom.Atom, atom2: Bio.PDB.Atom.Atom, atom3: Bio.PDB.Atom.Atom) -> HBondStats:
    assert _linalg_norm(atom0.coord - atom1.coord) <= 1.6, f"atoms 0,1 not bonded: {atom0.full_id} {atom1.full_id} ({np.linalg.norm(atom0.coord - atom1.coord)} > 1.6, {atom0.coord} {atom1.coord})"
    assert _linalg_norm(atom2.coord - atom3.coord) <= 1.6, f"atoms 2,3 not bonded: {atom2.full_id} {atom3.full_id} ({np.linalg.norm(atom0.coord - atom1.coord)} > 1.6, {atom0.coord} {atom1.coord})"
    assert _linalg_norm(atom1.coord - atom2.coord) > 1.3, f"atoms too close for h-bond: {atom1.full_id} {atom2.full_id} ({np.linalg.norm(atom1.coord - atom2.coord)} < 2, {atom1.coord} {atom2.coord})"
    assert _linalg_norm(atom1.coord - atom2.coord) < 10, f"atoms too far for h-bond: ({np.linalg.norm(atom1.coord - atom2.coord)} > 6, {atom1.coord} {atom2.coord})"
    return HBondStats(
        length=get_distance(atom1, atom2),
        donor_angle=get_angle(atom0, atom1, atom2),
        acceptor_angle=get_angle(atom3, atom2, atom1)
    )

resname_map = {
    'DT': 'U',
    'DC': 'C',
    'DA': 'A',
    'DG': 'G',
    'DU': 'U',
    'T': 'U',
}

def get_hbond_stats(pair_type: pair_defs.PairType, r1: AltResidue, r2: AltResidue) -> Optional[List[Optional[HBondStats]]]:
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
        hbond_stats(*atoms) if None not in atoms else None
        for atoms in bonds_atoms
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
                r1 = structure[model-1][chain1][(' ', res1, ins1 or ' ')]
                r1 = AltResidue(r1, alt1)
            except KeyError:
                print(f"Could not find residue1 {chain1}.{str(res1)+ins1} in {pdbid}")
                continue
            try:
                r2 = structure[model-1][chain2][(' ', res2, ins2 or ' ')]
                r2 = AltResidue(r2, alt2)
            except KeyError:
                print(f"Could not find residue2 {chain2}.{str(res2)+ins2} in {pdbid}")
                continue
            pair_type = pair_defs.PairType.create(pair_family, r1.resname, r2.resname, name_map=resname_map)
            hbonds = get_hbond_stats(pair_type, r1, r2)
            stats = None
            stats = calc_pair_stats(pair_type, r1, r2)
            if hbonds is not None:
                yield i, hbonds, stats
        except AssertionError as e:
            print(f"{e} on row:\n{to_csv_row(df_, i, 15)}")
            continue
        except Exception as e:
            print(f"Error on row:\n{to_csv_row(df_, i, 15)}")
            import traceback
            print(e)
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

def export_stats_csv(pdbid, df: pl.DataFrame, add_metadata_columns: bool, pair_type: Optional[str], max_bond_count: int) -> Tuple[str, pl.DataFrame, np.ndarray]:
    bond_params = [ x.name for x in dataclasses.fields(HBondStats) ]
    valid = np.zeros(len(df), dtype=np.bool_)
    columns: list[list[Optional[float]]] = [ [ None ] * len(df) for _ in range(max_bond_count * len(bond_params)) ]
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

    result_cols = {
        f"hb_{i}_{p}": pl.Series(c, dtype=pl.Float64)
        for i in range(max_bond_count)
        for p, c in zip(bond_params, columns[i * len(bond_params):])
    }
    result_cols.update(
        pl.DataFrame([
            p or PairStats(None, None, None, None, None)
            for p in pair_stats
        ]).to_dict()
    )
    result_df = pl.DataFrame(result_cols)
    if add_metadata_columns:
        h = structure.header if structure is not None else dict()
        result_df = result_df.with_columns(
            pl.lit(h.get('deposition_date', None), dtype=pl.Utf8).alias("deposition_date"),
            pl.lit(h.get('name', None), dtype=pl.Utf8).alias("structure_name"),
            pl.lit(h.get('structure_method', None), dtype=pl.Utf8).alias("structure_method"),
            pl.lit(h.get('resolution', None), dtype=pl.Float32).alias("resolution")
        )
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

def main(pool: Union[Pool, MockPool], args):
    df = load_inputs(pool, args)
    max_bond_count = get_max_bond_count(df)
    groups = list(df.groupby(pl.col("pdbid")))

    processes = []

    processes.append([
        pool.apply_async(export_stats_csv, args=[pdbid, group, args.metadata, args.pair_type, max_bond_count])
        for (pdbid,), group in groups
        # for chunk in group.iter_slices(n_rows=100)
    ])
    if args.dssr_binary is not None:
        import dssr_wrapper
        processes.append([
            pool.apply_async(dssr_wrapper.add_dssr_info, args=[pdbid, group, args.dssr_binary])
            for ((pdbid,), group) in groups
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

    df = pl.concat(result_chunks)
    if args.dedupe:
        df = remove_duplicate_pairs(df)
    df = df.sort('pdbid', 'model', 'nr1', 'nr2')
    df.write_csv(args.output)
    if not args.output.startswith("/dev/"):
        df.write_parquet(args.output + ".parquet")
    return df

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="""
        Adds geometric information to the specified pairing CSV files.
        Added columns:
            * hb_0_length, hb_0_heavy_a_angle, hb_1_length, hb_1_heavy_a_angle, hb_2_length, hb_2_heavy_a_angle - hydrogen bond lengths and angles between heavy atoms
            * coplanarity_angle - angle between planes of the two nucleotides
        When --dssr-binary is specified, DSSR --analyze is executed to gain additional information:
            * dssr_pairing_type - pairing type according to DSSR (e.g. WC, Platform, ~rHoogsteen)
            * pairing_type, shear, stretch, stagger, buckle, propeller, opening, shift, slide, rise, tilt, roll, twist
        """,
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("inputs", nargs="+")
    parser.add_argument("--pdbcache", nargs="+", help="Directories to search for PDB files in order to avoid downloading. Last directory will be written to, if the structure is not found and has to be downloaded.")
    parser.add_argument("--output", "-o", required=True, help="Output CSV file name")
    parser.add_argument("--threads", type=int, default=1, help="Maximum parallelism - number of worker processes to spawn")
    parser.add_argument("--metadata", type=bool, default=True, help="Add deposition_date, resolution and structure_method columns")
    parser.add_argument("--dssr-binary", type=str, help="If specified, DSSR --analyze will be invoked for each structure and its results stored as 'dssr_' prefixed columns")
    parser.add_argument("--filter", default=False, action="store_true", help="Filter out rows for which the values could not be calculated")
    parser.add_argument("--dedupe", default=False, action="store_true", help="Remove duplicate pairs, keep the one with shorter bonds or lower chain1,nr1")
    parser.add_argument("--pair-type", type=str, required=False)
    args = parser.parse_args()

    for x in args.pdbcache:
        pdb_utils.pdb_cache_dirs.append(os.path.abspath(x))
    os.environ["PDB_CACHE_DIR"] = ';'.join(pdb_utils.pdb_cache_dirs)

    multiprocessing.set_start_method("spawn")
    if args.threads == 1:
        main(MockPool(), args)
    else:
        with multiprocessing.Pool(processes=args.threads) as pool:
            main(pool, args)
