import polars as pl
import argparse, os
import pair_defs
import pair_csv_parse
import simulate_fr3d_output

def write_directory(df, out, v):
    if {'pdbid', 'model', 'chain1', 'chain2', 'nr1', 'nr2', 'family'}.issubset(df.columns):
        if v: print("basepair.txt...", end="")
        identification_df = df.select(
                pl.col("family"),
                pl.col("pdbid"),
                pl.col("model"),
                pl.col(r"^(chain|nr|res|alt|ins|symmetry_operation)\d$")
            ).collect()
        simulate_fr3d_output.write_file(
            os.path.join(out, "basepair.txt"),
            pdbid_global=None,
            df=identification_df,
            detailed=True,
            only_once=False,
            comparison_column=False,
            additional_columns=[]
        )
        if v: print(" Done.")
        if v: print("id.csv...", end="")

        identification_df.write_csv(os.path.join(out, "id.csv"), include_header=True)
        if v: print(" Done.")

        del identification_df

    numeric_columns = [c for c, type in df.schema.items() if type.is_float()]

    for c in numeric_columns:
        if c in ["mode_deviations"]:
            continue
        if v: print(f"{c}.csv...", end="")
        df.select(pl.col(c)).sink_csv(os.path.join(out, f"{c}.csv"), include_header=False)
        if v: print(f" Done.")

def main(args):
    df = pair_csv_parse.scan_pair_csvs(args.input)
    if "res1" in df.columns:
        df = df.with_columns(
            res1_=pl.col("res1").replace(pair_defs.resname_map).str.to_uppercase(),
            res2_=pl.col("res2").replace(pair_defs.resname_map).str.to_uppercase()
        )
    out = args.output_dir
    if not args.quiet:
        print("schema:")
        print(', '.join(f"{c}:{t}"  for c, t in df.schema.items()))
    os.makedirs(out, exist_ok=True)
    if os.listdir(out):
        print(f"WARNING: Output directory {args.output_dir} is not empty")

    if not args.partition:
        write_directory(df, out, not args.quiet)
    else:
        for key, g in df.collect().group_by(args.partition.split(",")):
            out_p = os.path.join(out, "-".join(key))
            os.makedirs(out_p, exist_ok=True)
            write_directory(g.lazy(), out_p, not args.quiet)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('input', nargs="+", type=str)
    parser.add_argument('output_dir', type=str)
    parser.add_argument('--partition', type=str)
    parser.add_argument('--quiet', '-q', default=False, action="store_true")
    args = parser.parse_args()
    main(args)
