import functools
import itertools
from typing import Any, List, Optional, Tuple
import Bio.PDB, Bio.PDB.Structure, Bio.PDB.Model, Bio.PDB.Residue, Bio.PDB.Atom, Bio.PDB.Chain, Bio.PDB.Entity
import os, sys, io, gzip, math, numpy as np
import polars as pl
from dataclasses import dataclass
import dataclasses
from requests import HTTPError
import multiprocessing
from pair_csv_parse import scan_pair_csvs
import pdb_utils

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
    if sum_error > 1 * len(planar_atoms): # TODO: reasonable threshold?
        raise ResideTransformError(f"Residue really doesn't seem planar: {res.res.full_id}, plane RMSE = {np.sqrt(sum_error / len(planar_atoms))}")

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
    bogopropeller: float
    # opening: float
    # nearest_distance: float

def calc_pair_stats(res1: AltResidue, res2: AltResidue) -> Optional[PairStats]:
    p1 = try_get_residue_posinfo(res1)
    p2 = try_get_residue_posinfo(res2)
    if p1 is None or p2 is None:
        return None
    res1_ = AltResidue(transform_residue(res1.res, p1.translation, p1.rotation), res1.alt)
    res2_ = AltResidue(transform_residue(res2.res, p1.translation, p1.rotation), res2.alt)

    bogopropeller = np.arccos(np.dot(p1.plane_normal_vector, p2.plane_normal_vector))
    return PairStats(bogopropeller)

def get_angle(atom1: Bio.PDB.Atom.Atom, atom2: Bio.PDB.Atom.Atom, atom3: Bio.PDB.Atom.Atom) -> float:
    v1 = atom1.coord - atom2.coord
    v2 = atom3.coord - atom2.coord
    return np.arccos(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2)))

def get_distance(atom1: Bio.PDB.Atom.Atom, atom2: Bio.PDB.Atom.Atom) -> float:
    return float(np.linalg.norm(atom1.coord - atom2.coord))


@dataclass
class HBondStats:
    length: float
    donor_angle: float
    acceptor_angle: float

def hbond_stats(atom0: Bio.PDB.Atom.Atom, atom1: Bio.PDB.Atom.Atom, atom2: Bio.PDB.Atom.Atom, atom3: Bio.PDB.Atom.Atom) -> HBondStats:
    assert np.linalg.norm(atom0.coord - atom1.coord) <= 1.6, f"atoms 0,1 not bonded: {atom0.full_id} {atom1.full_id} ({np.linalg.norm(atom0.coord - atom1.coord)} > 1.6, {atom0.coord} {atom1.coord})"
    assert np.linalg.norm(atom2.coord - atom3.coord) <= 1.6, f"atoms 2,3 not bonded: {atom2.full_id} {atom3.full_id} ({np.linalg.norm(atom0.coord - atom1.coord)} > 1.6, {atom0.coord} {atom1.coord})"
    assert np.linalg.norm(atom1.coord - atom2.coord) > 1.3, f"atoms too close for h-bond: {atom1.full_id} {atom2.full_id} ({np.linalg.norm(atom1.coord - atom2.coord)} < 2, {atom1.coord} {atom2.coord})"
    assert np.linalg.norm(atom1.coord - atom2.coord) < 6, f"atoms too far for h-bond: ({np.linalg.norm(atom1.coord - atom2.coord)} > 6, {atom1.coord} {atom2.coord})"
    return HBondStats(
        length=get_distance(atom1, atom2),
        donor_angle=get_angle(atom0, atom1, atom2),
        acceptor_angle=get_angle(atom3, atom2, atom1)
    )

def hbond_swap_nucleotides(hbond: tuple[str, str, str, str]) -> tuple[str, str, str, str]:
    return tuple(
        atom.replace("A", "ⒷⒷⒷⒷⒷⒷ").replace("B", "A").replace("ⒷⒷⒷⒷⒷⒷ", "B")
        for atom in hbond
    ) # type: ignore

hbonding_atoms: dict[tuple[str, str], list[tuple[str, str, str, str]]] = {
    ('cWW', 'C-G'): [
        ('AC4', 'AN4', 'BO6', 'BC6'),
        ('BC6', 'BN1', 'AN3', 'AC2'),
        ('BC2', 'BN2', 'AO2', 'AC2'),
    ],
    ('cWW', 'A-U'): [
        ('AC6', 'AN6', 'BO4', 'BC4'),
        ('BC2', 'BN3', 'AN1', 'AC2')
    ],
    ('cWW', 'G-U'): [
        ('BC4', 'BN3', 'AO6', 'AC6'),
        ('AC2', 'AN1', 'BO2', "BC2"),
        ('AC2', 'AN2', "BO2'", "BC2'"), # ??? pres vodu
    ],
    ('tHS', 'A-G'): [
        ('BC2', 'BN2', 'AN7', 'AC8'),
        ('AC6', 'AN6', 'BN3', 'BC4'),
        ('AC6', 'AN6', "BO2'", "BC2'"),
    ],
    # ('cSS', 'A-C'): [ ], # kokotina
    ('tSS', 'A-G'): [
        ("BC2", "BN2", "AN3", "AC2"),
        ("AN3", "AC2", "BN3", "BC4"),
        ("BC2'", "BO2'", "AN1", "AC6")
    ],
    ('tHW', 'A-U'): [
        ("BC4", "BN3", "AN7", "AC8"),
        ("AC6", "AN6", "BO2", "BC2"),
    ],
    ('cWW', 'U-U'): [
        ("BC4", "BN3", "AO4", "AC4"),
        ("AC4", "AN3", "BO2", "BC2"),
        ("BC2'", "BO2'", "AO2", "AC2"), # pres vodu
    ],
    ('tHH', 'A-A'): [
        ("AC6", "AN6", "BN7", "BC8"),
        ("BC6", "BN6", "AN7", "AC8"),
        ("BC6", "BN6", "AOP2", "AP"),
    ],
    ('cSS', 'A-G'): [
        ("BC2", "BN2", "AN3", "AC4"),
        ("AC2'", "AO2'", "BN3", "BC4"),
        ("BC2'", "BO2'", "AO2'", "AC2'"),
    ],
    ('cWW', 'A-G'): [
        ('AC6', 'AN6', 'BO6', 'BC6'),
        ('BC6', 'BN1', 'AN1', 'AC6'),
        # ('BC2', 'BN2', 'AC2', 'AN3'),
    ],
    ('tWW', 'A-A'): [
        ("AC6", "AN6", "BN1", "BC6"),
        ("BC6", "BN6", "AN1", "AC6")
    ],
    # ('tSS', 'A-C'): [ ], # kokotina
    ('tWS', 'A-G'): [
        ("BC2", "BN2", "AN1", "AC6"),
        ("AC6", "AN6", "BN3", "BC4"),
        ("AC6", "AN6", "BO2'", "BC2'"),
    ],
    # ('cSS', 'A-U'): [ ], # kokotina
    # ('cSH', 'G-U'): [ ], # kokotina
    ('cSS', 'A-A'): [
        ("BN3", "BC2", "AN3", "AC4"),
        ("AC2'", "AO2'", "BN3", "BC4"),
        ("BC2'", "BO2'", "AO2'", "AC2'"),
    ],
    ('tSS', 'G-G'): [
        ("BC2", "BN2", "AN3", "AC4"),
        ("AC2", "AN2", "BN3", "BC4"),
        ("AC2", "AN2", "BO2'", "BC2'"),
        # ("BC2", "BN2", "AO2'", "AC2'"),
    ],
    ('cWH', 'G-G'): [
        ("AC2", "AN2", "BN7", "BC8"),
        ("AC6", "AN1", "BO6", "BC6")
    ],
    ('cWW', 'C-U'): [
        ("AC4", "AN4", "BO4", "BC4"),
        ("BC4", "BN3", "AN3", "AC4"),
    ],
    # ('cWH', 'C-U'): [
    #     ("AC4", "AN4", "BO4", "BC4"),
    #     ("BC2", "BN3", "AN3", "AC4"),
    # ],
}

resname_map = {
    'DT': 'U',
    'DC': 'C',
    'DA': 'A',
    'DG': 'G',
    'DU': 'U',
    'T': 'U',
}

def get_hbond_stats(pair_type: str, r1: AltResidue, r2: AltResidue) -> Optional[List[Optional[HBondStats]]]:
    if r1.resname > r2.resname:
        r1, r2 = r2, r1
    pair_name = resname_map.get(r1.resname, r1.resname) + '-' + resname_map.get(r2.resname, r2.resname)
    if (pair_type, pair_name) not in hbonding_atoms:
        return None
    def get_atom(n):
        if n[0] == 'A':
            return r1.get_atom(n[1:], None)
        elif n[0] == 'B':
            return r2.get_atom(n[1:], None)
        else:
            raise ValueError(f"Invalid atom name: {n}")

    bonds_atoms = [
        tuple(get_atom(a) for a in atoms)
        for atoms in hbonding_atoms[(pair_type, pair_name)]
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


def get_stats_for_csv(df_: pl.DataFrame, structure: Bio.PDB.Structure.Structure, pair_type: str):
    df = df_.with_row_count()
    for (i, pdbid, model, chain1, res1, ins1, alt1, chain2, res2, ins2, alt2) in zip(df["row_nr"], df["pdbid"], df["model"], df["chain1"], df["nr1"], df["ins1"], df["alt1"], df["chain2"], df["nr2"], df["ins2"], df["alt2"]):
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
            hbonds = get_hbond_stats(pair_type, r1, r2)
            # stats = None
            stats = calc_pair_stats(r1, r2)
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
    def pair_id(chain1, res1, ins1, alt1, chain2, res2, ins2, alt2):
        pair = [ (chain1, res1, ins1, alt1), (chain2, res2, ins2, alt2) ]
        pair.sort()
        return "|".join(tuple(str(x) for x in (pair[0] + pair[1])))
    pair_ids = [ pair_id(chain1, res1, ins1, alt1, chain2, res2, ins2, alt2) for chain1, res1, ins1, alt1, chain2, res2, ins2, alt2 in zip(df["chain1"], df["nr1"], df["ins1"], df["alt1"], df["chain2"], df["nr2"], df["ins2"], df["alt2"]) ]
    df = df.with_columns(
        pl.Series(pair_ids, dtype=pl.Utf8).alias("_tmp_pair_id")
    ).with_row_count("_tmp_row_nr")
    score = pl.lit(0, dtype=pl.Float64)
    for col in df.columns:
        if col.startswith("hb_") and col.endswith("_length"):
            score += pl.col(col).fill_null(100)
        if col.startswith("dssr_"):
            score += pl.col(col).is_null().cast(pl.Float64) * 3
        if col == "bogopropeller":
            score += pl.col(col).is_null().cast(pl.Float64) * 10

    df = df.sort([score, "chain1", "nr1", "ins1", "alt1", "chain2", "nr2", "ins2", "alt2"])
    df = df.unique(["pdbid", "model", "_tmp_pair_id"], keep="first", maintain_order=True)
    df = df.sort("_tmp_row_nr")
    df = df.drop(["_tmp_pair_id", "_tmp_row_nr"])
    return df

def export_stats_csv(pdbid, df: pl.DataFrame, add_metadata_columns: bool, pair_type: str) -> Tuple[str, pl.DataFrame, np.ndarray]:
    bond_count = 3
    bond_params = [ x.name for x in dataclasses.fields(HBondStats) ]
    valid = np.zeros(len(df), dtype=np.bool_)
    columns: list[list[Optional[float]]] = [ [ None ] * len(df) for _ in range(bond_count * len(bond_params)) ]
    bogopropeller: list[Optional[float]] = [ None ] * len(df)

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

            if stats:
                bogopropeller[i] = stats.bogopropeller

    result_cols = {
        f"hb_{i}_{p}": pl.Series(c, dtype=pl.Float64) for i in range(bond_count) for p, c in zip(bond_params, columns[i * len(bond_params):])
    }
    result_cols["bogopropeller"] = pl.Series(bogopropeller, dtype=pl.Float64)
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

def main(pool, args):
    df = scan_pair_csvs(args.csvs).sort('pdbid', 'model', 'nr1', 'nr2').collect()

    groups = list(df.groupby(pl.col("pdbid")))

    processes = []

    processes.append([
        pool.apply_async(export_stats_csv, args=[pdbid, group, args.metadata, args.pair_type])
        for pdbid, group in groups
        # for chunk in group.iter_slices(n_rows=100)
    ])
    if args.dssr_binary is not None:
        import dssr_wrapper
        processes.append([
            pool.apply_async(dssr_wrapper.add_dssr_info, args=[pdbid, group, args.dssr_binary])
            for (pdbid, group) in groups
        ])

    result_chunks = []
    for (_pdbid, group), p in zip(groups, zip(*processes)):
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
    df.write_parquet(args.output + ".parquet")
    return df

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="""
        Adds geometric information to the specified pairing CSV files.
        Added columns:
            * hb_0_length, hb_0_heavy_a_angle, hb_1_length, hb_1_heavy_a_angle, hb_2_length, hb_2_heavy_a_angle - hydrogen bond lengths and angles between heavy atoms
            * bogopropeller - angle between planes of the two nucleotides
        When --dssr-binary is specified, DSSR --analyze is executed to gain additional information:
            * dssr_pairing_type - pairing type according to DSSR (e.g. WC, Platform, ~rHoogsteen)
            * pairing_type, shear, stretch, stagger, buckle, propeller, opening, shift, slide, rise, tilt, roll, twist
        """,
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("csvs", nargs="+")
    parser.add_argument("--pdbcache", nargs="+", help="Directories to search for PDB files in order to avoid downloading. Last directory will be written to, if the structure is not found and has to be downloaded.")
    parser.add_argument("--output", "-o", required=True, help="Output CSV file name")
    parser.add_argument("--threads", type=int, default=1, help="Maximum parallelism - number of worker processes to spawn")
    parser.add_argument("--metadata", type=bool, default=True, help="Add deposition_date, resolution and structure_method columns")
    parser.add_argument("--dssr-binary", type=str, help="If specified, DSSR --analyze will be invoked for each structure and its results stored as 'dssr_' prefixed columns")
    parser.add_argument("--filter", default=False, action="store_true", help="Filter out rows for which the values could not be calculated")
    parser.add_argument("--dedupe", default=False, action="store_true", help="Remove duplicate pairs, keep the one with shorter bonds or lower chain1,nr1")
    parser.add_argument("--pair-type", default="cWW")
    args = parser.parse_args()

    for x in args.pdbcache:
        pdb_utils.pdb_cache_dirs.append(os.path.abspath(x))
    os.environ["PDB_CACHE_DIR"] = ';'.join(pdb_utils.pdb_cache_dirs)

    multiprocessing.set_start_method("spawn")
    with multiprocessing.Pool(processes=args.threads) as pool:
        main(pool, args)
