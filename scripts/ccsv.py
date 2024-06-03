import polars as pl
import argparse, os
import pair_csv_parse
import simulate_fr3d_output

def main(args):
    df = pair_csv_parse.scan_pair_csvs(args.input)
    out = args.output_dir
    print("schema:")
    print(', '.join(f"{c}:{t}"  for c, t in df.schema.items()))
    os.makedirs(out, exist_ok=True)
    if os.listdir(out):
        print(f"WARNING: Output directory {args.output_dir} is not empty")

    if {'pdbid', 'model', 'chain1', 'chain2', 'nr1', 'nr2', 'family'}.issubset(df.columns):
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
        print("Done: basepair.txt")

        identification_df.write_csv(os.path.join(out, "id.csv"), include_header=True)
        print("Done: id.csv")

        del identification_df

    numeric_columns = [c for c, type in df.schema.items() if type.is_float()]

    for c in numeric_columns:
        df.select(pl.col(c)).sink_csv(os.path.join(out, f"{c}.csv"), include_header=False, null_value="NaN")
        print(f"Done: {c}.csv")

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('input', nargs="+", type=str)
    parser.add_argument('output_dir', type=str)
    args = parser.parse_args()
    main(args)
