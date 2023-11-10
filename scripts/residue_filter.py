import dataclasses
from typing import Optional
import polars as pl
import numpy as np
import os, re, math


@dataclasses.dataclass
class NucleotideID:
    pdbid: str
    chain: str
    base: str
    id: int
    ins: Optional[str]
    alt: Optional[str]
def read_id_list(file) -> pl.DataFrame:
    results = []
    with open(file) as f:
        for l in f:
            m = re.match(r"^(?P<pdbid>[a-zA-Z0-9]{4})_(?P<chain>[a-zA-Z0-9]+)_(?P<base>[a-zA-Z0-9]+)(\.(?P<alt>[a-zA-Z0-9]+))?_(?P<id>-?\d+)(\.(?P<ins>[a-zA-Z0-9]+))?$", l)
            if m is None:
                raise Exception(f"Could not parse id: {l}")
            results.append(NucleotideID(
                pdbid=m.group("pdbid"),
                chain=m.group("chain"),
                base=m.group("base"),
                id=int(m.group("id")),
                ins=m.group("ins") or None,
                alt=m.group("alt") or None,
            ))
    return pl.DataFrame(results, schema={ "pdbid": pl.Utf8, "chain": pl.Utf8, "base": pl.Utf8, "id": pl.Int64, "ins": pl.Utf8, "alt": pl.Utf8 })

def read_id_lists(directory) -> dict[str, pl.DataFrame]:
    return {
        "DNA-0-1.8": read_id_list(directory + "/non_redundant_DNA_filtered_residues_0-1.8A.csv").with_columns(pl.lit(True).alias("dna_hres")),
        "RNA-0-1.8": read_id_list(directory + "/non_redundant_RNA_filtered_residues_0-1.8A.csv").with_columns(pl.lit(True).alias("rna_hres")),
        "DNA-1.8-3.5": read_id_list(directory + "/non_redundant_DNA_filtered_residues_1.8-3.5A.csv").with_columns(pl.lit(True).alias("dna_lres")),
        "RNA-1.8-3.5": read_id_list(directory + "/non_redundant_RNA_filtered_residues_1.8-3.5A.csv").with_columns(pl.lit(True).alias("rna_lres")),
    }

def add_res_filter_columns(df, residue_lists: dict[str, pl.DataFrame]):
    print(next(iter(residue_lists.values())))
    rcols = { f"{name}-r{resix}":
        df.join(reslist.with_columns(pl.lit(True).alias("__tmp")), left_on=["pdbid", f"chain{resix}", f"res{resix}", f"nr{resix}", f"alt{resix}", f"ins{resix}"], right_on=["pdbid", "chain", "base", "id", "alt", "ins"], how="left")
            .get_column("__tmp")
            .fill_null(False)
        for name, reslist in residue_lists.items()
        for resix in [1, 2]
    }
    df = df.with_columns(**rcols)
    df = df.with_columns(**{
        name: df[f"{name}-r1"] | df[f"{name}-r2"]
        for name in residue_lists.keys()
    })
    label_condition = pl
    for name in residue_lists.keys():
        label_condition = label_condition.when(pl.col(name)).then(pl.lit(name))
    label_condition = \
        label_condition.when(pl.col("res1").str.starts_with("D")).then(pl.lit("DNA")) \
        .otherwise(pl.lit("RNA"))
    df = df.with_columns(label_condition.alias("label"))
    return df

def _read_inputs(files):
    import pair_csv_parse
    return pl.read_parquet(files[0]) if files[0].endswith(".parquet") else pair_csv_parse.scan_pair_csvs(files, header=False).collect()

def _stats(df, residue_lists: dict[str, pl.DataFrame]):
    results = {}
    for name in residue_lists.keys():
        results[name] = df.get_column(name).sum()

    r1dna = pl.col("res1").str.starts_with("D")
    r2dna = pl.col("res2").str.starts_with("D")
    results["all_RNA"] = len(df.filter(~r1dna & ~r2dna))
    results["all_DNA: "] = len(df.filter(r1dna & r2dna))
    results["all_mixed"] = len(df.filter((r1dna & ~r2dna) | (~r1dna & r2dna)))
    results["total"] = len(df)
    return results

def _run_filter(args):
    df = _read_inputs(args.input)
    residue_lists = read_id_lists(args.residue_directory)
    df = add_res_filter_columns(df, residue_lists)
    if args.output is not None:
        if args.output.endswith(".parquet"):
            df.write_parquet(args.output)
        else:
            df.write_csv(args.output)

    if args.stats:
        s = _stats(df, residue_lists)
        for k, v in s.items():
            print(f"{k}: {v}")

def _run_stats(args):
    files = os.listdir(args.input_directory)
    pairing_types = [ f for f in files if re.match(r"^[ct][HSW]{2}", f) ]
    if len(pairing_types) == 0:
        pairing_types = ["."]

    residue_lists = read_id_lists(args.residue_directory)
    stats = []
    for pairing_type in pairing_types:
        pdir = os.path.join(args.input_directory, pairing_type)
        print(pdir)
        for pairing in [("D?A", "D?A"), ("D?A", "D?C"), ("D?A", "D?G"), ("D?A", "D?[TU]"), ("D?C", "D?C"), ("D?C", "D?G"), ("D?C", "D?[TU]"), ("D?G", "D?G"), ("D?G", "D?[TU]"), ("D?[TU]", "D?[TU]")]:
            files = [ os.path.join(pdir, f) for f in os.listdir(pdir) if re.match(f"^{pairing[0]}[_-]{pairing[1]}\\.(csv|parquet)$", f) ]
            files = [ f for f in files if os.path.isfile(f) and os.stat(f).st_size > 0 ]
            print(files)
            if len(files) == 0:
                continue

            df = _read_inputs(files)
            df = add_res_filter_columns(df, residue_lists)
            s = _stats(df, residue_lists)
            stats.append({
                "pairing_type": pairing_type,
                "pairing": f"{pairing[0]}-{pairing[1]}",
                **s
            })
    sdf = pl.DataFrame(stats)
    sdf = sdf.sort(pl.col("DNA-0-1.8") + pl.col("RNA-0-1.8"), descending=True)
    if args.output is not None:
        if args.output.endswith(".parquet"):
            sdf.write_parquet(args.output)
        else:
            sdf.write_csv(args.output)

    else:
        pl.Config.set_tbl_cols(100)
        pl.Config.set_tbl_rows(1000000)
        print(sdf)
    with open("tablica.html", "wt") as f:
        f.write(sdf._repr_html_())

    

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--residue_directory", required=True)
    subs = parser.add_subparsers(required=True)
    filter_cmd = subs.add_parser('filter')
    filter_cmd.add_argument("input", nargs="+")
    filter_cmd.add_argument("--output", default=None)
    filter_cmd.add_argument("--stats", action="store_true")
    filter_cmd.set_defaults(func=_run_filter)
    stats_cmd = subs.add_parser('stats')
    stats_cmd.add_argument("input_directory")
    stats_cmd.add_argument("--output", default=None)
    stats_cmd.set_defaults(func=_run_stats)
    args = parser.parse_args()
    args.func(args)



