import multiprocessing
import subprocess, os, sys, re, math
from typing import List, Literal, Optional, Sequence, TextIO
import numpy as np
import polars as pl
import dataclasses
from dataclasses import dataclass
import tempfile
import pdb_utils

@dataclass(frozen=True, order=True)
class NucleotideID:
    sequential_id: int
    chain: str
    id: str
    base: str
    alt: str
    def replace(self, **kwargs):
        return dataclasses.replace(self, **kwargs)

@dataclass(frozen=True, order=True)
class Pair:
    nt1: NucleotideID
    nt2: NucleotideID
    pairing_type: str
    shear: Optional[float]
    stretch: Optional[float]
    stagger: Optional[float]
    buckle: Optional[float]
    propeller: Optional[float]
    opening: Optional[float]
    shift: Optional[float]
    slide: Optional[float]
    rise: Optional[float]
    tilt: Optional[float]
    roll: Optional[float]
    twist: Optional[float]
    def replace(self, **kwargs):
        return dataclasses.replace(self, **kwargs)

@dataclass
class AnalysisResult:
    pairs: list[Pair]

def _skip_comments(lines: Sequence[str]):
    return (line for line in lines if not line.startswith('#'))

def _parse_id(x: str) -> tuple[str, str, str, str]:
    m = re.match(r"^(?P<chain>\w+)\.(?P<base>\w+?)/?(?P<id>-?\d+)(\^(?P<alt>\w+))?$", x)
    if m is None:
        raise Exception(f"Could not parse DSSR id: {x}")
    
    return m.group("chain"), m.group("id"), m.group("base"), m.group("alt") or ''

def _normalize_base(b):
    b = b.upper()
    if len(b) > 1 and b.startswith('D'):
        b = b[1:]
    if b == '2MG' or b == '7MG' or b == 'OMG' or b == 'GTP':
        return 'G'
    if b == '5MC' or b == 'OMC':
        return 'C'
    if b == '5MU' or b == 'CTG':
        return 'T'
    if b == 'A2M':
        return 'A'

    return b

def _read_output_files(directory: str) -> list[Pair]:
    with open(os.path.join(directory, 'dssr-pairs.txt')) as f:
        pairs = list(_skip_comments(f.readlines()))
    with open(os.path.join(directory, 'dssr-dsStepPars.txt')) as f:
        parameters = list(_skip_comments(f.readlines()))

    result = []
    for pair, pars in zip(pairs, parameters):
        seq1, seq2, wtf, wtf, index_asi, id1, id2, base_pair, pairing_type = pair.split()
        base1_, base2_ = re.split(r'[-+]', base_pair)
        chain1, id1, base1, alt1 = _parse_id(id1)
        chain2, id2, base2, alt2 = _parse_id(id2)
        # assert base1 == base1_, f"{base1} != {base1_} in {pair} / {pars}"
        # assert base2 == base2_, f"{base2} != {base2_} in {pair} / {pars}"
        if len(_normalize_base(base1)) == 1 and len(_normalize_base(base2)) == 1:
            # ^ don't spam about modified bases 
            if _normalize_base(base1) != base1_.upper() or _normalize_base(base2) != base2_.upper():
                print(f"WARNING: base mismatch in ({base1}/{base1_}, {base2}/{base2_}) {pair} / {pars}")
                continue

        base_pair_, *parameters = pars.split()
        assert base_pair == base_pair_
        parameters = [float(x) if x != "999999" else None for x in parameters]

        result.append(Pair(
            NucleotideID(int(seq1), chain1, id1, base1, alt1),
            NucleotideID(int(seq2), chain2, id2, base2, alt2),
            pairing_type,
            *parameters
        ))
    return result

def run_dssr(dssr_binary, pdbid, input_file: TextIO) -> list[Pair]:
    with tempfile.TemporaryDirectory(f'-{pdbid}-dssr') as dir:
        with open(os.path.join(dir, 'input.cif'), 'w') as f:
            with input_file as f2:
                f.write(f2.read())
        try:
            subprocess.run([dssr_binary, '--analyze', '-i=' + os.path.join(dir, 'input.cif')], check=True, cwd=dir,
                           stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL
                        )
        except subprocess.CalledProcessError as e:
            print(f"WARNING: dssr failed for {pdbid}: {e}")
            return []
        if not os.path.exists(os.path.join(dir, 'dssr-pairs.txt')):
            print(f"WARNING: no dssr output for {pdbid}")
            return []
        try:
            return _read_output_files(dir)
        except Exception as e:
            print(f"WARNING: dssr output parsing failed for {pdbid}: {e}")
            import traceback
            traceback.print_exc()
            return []

def add_dssr_info(pdbid, df: pl.DataFrame, dssr_binary, column_prefix = 'dssr_') -> tuple[str, pl.DataFrame, np.ndarray]:
    pairs = run_dssr(dssr_binary, pdbid, pdb_utils.open_pdb_file(None, pdbid))
    lookup = { (pair.nt1.replace(sequential_id=-1), pair.nt2.replace(sequential_id=-1)): pair for pair in pairs }
    # for p in pairs:
    #     print(p)
    def lookup_pair(nt1: NucleotideID, nt2: NucleotideID) -> Pair | None:
        default = None
        if nt1.alt or nt2.alt:
            default = lookup_pair(nt1.replace(alt=''), nt2.replace(alt=''))
        return lookup.get((nt1, nt2), lookup.get((nt2, nt1), default))

    columns_names = [ "pairing_type", "shear", "stretch", "stagger", "buckle", "propeller", "opening", "shift", "slide", "rise", "tilt", "roll", "twist" ]
    columns = [ [None] * len(df) for _ in columns_names ]
    for i, (pdbid_, model, chain1, res1, nr1, ins1, alt1, chain2, res2, nr2, ins2, alt2) in enumerate(df[['pdbid', 'model', 'chain1', 'res1', 'nr1', 'ins1', 'alt1', 'chain2', 'res2', 'nr2', 'ins2', 'alt2']].iter_rows()):
        if pdbid_ != pdbid:
            print(f"WARNING: pdbid mismatch: {pdbid_} != {pdbid}")

        if alt1 == '?':
            alt1 = None
        if alt2 == '?':
            alt2 = None
        if ins1.strip() == '':
            ins1 = None
        if ins2.strip() == '':
            ins2 = None

        pair = lookup_pair(
            NucleotideID(-1, chain1, str(nr1), res1, alt1 or ins1 or ''),
            NucleotideID(-1, chain2, str(nr2), res2, alt2 or ins2 or ''))
        if pair is None:
            print(f"WARNING: no pair for {pdbid} {chain1}_{res1}_{nr1}_{ins1}_{alt1} {chain2}_{res2}_{nr2}_{ins2}_{alt2}")
            continue

        for col, colname in zip(columns, columns_names):
            col[i] = getattr(pair, colname)

    assert len(columns[0]) == len(df), f"{len(columns[0])} != {len(df)} for {pdbid}, col0 = {columns[0]}"
    result_columns = {
        column_prefix + name: pl.Series(column, dtype=pl.Utf8 if name == "pairing_type" else pl.Float64)
        for column, name in zip(columns, columns_names)
    }
    valid = np.array([ x != None for x in columns[0] ], dtype=np.bool_)
    if np.sum(valid) != len(df):
        print(f"WARNING: {len(df) - np.sum(valid)} pairs out of {len(df)} in {pdbid} don't match any DSSR pair. Total DSSR pairs = {len(pairs)}")
    return pdbid, pl.DataFrame(result_columns), valid
    

if __name__ == "__main__":
    import argparse, pairs
    parser = argparse.ArgumentParser()
    parser.add_argument("csvs", nargs="+")
    parser.add_argument("--pdbcache", nargs="+")
    parser.add_argument("--output", "-o", required=True)
    parser.add_argument("--dssr-binary", required=True)
    parser.add_argument("--threads", type=int, default=1)
    args = parser.parse_args()

    for x in args.pdbcache:
        pdb_utils.pdb_cache_dirs.append(os.path.abspath(x))
    os.environ["PDB_CACHE_DIR"] = ';'.join(pdb_utils.pdb_cache_dirs)

    multiprocessing.set_start_method("spawn")
    with multiprocessing.Pool(processes=args.threads) as pool:
        df = pairs.scan_pair_csvs(args.csvs).sort('pdbid', 'model', 'nr1', 'nr2').collect()

        procs = [
            pool.apply_async(add_dssr_info, args=[pdbid, group, args.dssr_binary])
            for pdbid, group in df.groupby(pl.col("pdbid"))
        ]
        df = pl.concat([ p for p in [ p.get()[0] for p in procs ] if len(p) > 0 ])
        df.sort('pdbid', 'model', 'nr1', 'nr2').write_csv(args.output)
