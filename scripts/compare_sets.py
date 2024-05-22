#!/usr/bin/env python3

import typing as ty
import os, sys, io, re, gzip, math, json, functools, itertools, numpy as np
import polars as pl
from dataclasses import dataclass
import dataclasses
from pair_csv_parse import scan_pair_csvs
import pair_defs as pair_defs
import duckdb
import subprocess


def main(args):
    boundaries = pl.read_csv(args.boundaries, infer_schema_length=100000) if args.boundaries else None

    if args.comparison:
        comparison_file = args.comparison
        file_is_tmp = False
    else:
        if args.output_full_diff:
            comparison_file = args.output_full_diff
            file_is_tmp = False
        else:
            comparison_file = "tmp-duckdb-union.parquet"
            file_is_tmp = True

        pairid_sql = "(pdbid || '-' || model || '-' || chain1 || '_' || coalesce(alt1, '') || nr1 || coalesce(ins1, '') || '-' || chain2 || '_' || coalesce(alt2, '') || nr2 || coalesce(ins2, ''))"
        duckdb.execute(f"""
            create view current AS SELECT {pairid_sql} as pairid, * FROM parquet_scan('{args.target}')
            """, [])
        duckdb.execute(f"""
            create view baseline AS SELECT {pairid_sql} as pairid, * FROM parquet_scan('{args.baseline}')
            """, [])
        duckdb.execute(f"""
            COPY (SELECT DISTINCT ON (pairid, family)
                * EXCLUDE (comparison_in_baseline, comparison_in_current),
                bool_or(comparison_in_baseline) OVER (PARTITION BY pairid) as comparison_in_baseline,
                bool_or(comparison_in_current) OVER (PARTITION BY pairid) as comparison_in_current
                FROM (
                    SELECT *, TRUE AS comparison_in_current, FALSE AS comparison_in_baseline FROM current
                    UNION ALL BY NAME
                    SELECT *, FALSE AS comparison_in_current, TRUE AS comparison_in_baseline FROM baseline
                    ORDER BY comparison_in_current DESC
                )) to '{comparison_file}' (format parquet, COMPRESSION 'zstd');
            """)
        
    
    comparison = pl.scan_parquet(comparison_file, cache=False, low_memory=True)
    comparison = comparison.with_columns(
        bases=pl.col("res1").replace(pair_defs.resname_map).str.to_uppercase() + "-" + pl.col("res2").replace(pair_defs.resname_map).str.to_uppercase(),
    )
    comparison = comparison.group_by([ pl.col("family").str.to_lowercase(), pl.col("bases") ]).agg(
        count_all=pl.len(),
        count_baseline=pl.col("comparison_in_baseline").fill_null(False).sum(),
        count_target=pl.col("comparison_in_current").fill_null(False).sum(),
        count_dropped=(pl.col("comparison_in_baseline").fill_null(False) & ~pl.col("comparison_in_current").fill_null(False)).sum(),
        count_added=(~pl.col("comparison_in_baseline").fill_null(False)  & pl.col("comparison_in_current").fill_null(False)).sum(),
    ).collect(streaming=True)

    comparison = comparison.with_columns(
        family_id=pl.col("family").replace(pair_defs.pair_families_ids, default=None, return_dtype=pl.Int32),
        diff = pl.format("-{} +{}", pl.col("count_dropped"), pl.col("count_added")),
        diff_percent = pl.format("-{}% +{}%", (pl.col("count_dropped") / pl.col("count_baseline") * 100).round_sig_figs(2), (pl.col("count_added") / pl.col("count_target") * 100).round_sig_figs(2)),
    )

    if boundaries is not None:
        comparison = boundaries.select("family", "bases").join(comparison, on=[pl.col("family").str.to_lowercase(), "bases"], how="left")
    else:
        comparison = comparison.sort("family_id", "bases")
    
    # comparison = comparison.select([
    #     "family", "bases", "count_baseline", "count_target", "count_dropped", "count_added", "diff",
    # ])
    comparison.write_csv(args.output)

    if file_is_tmp:
        os.remove(comparison_file)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="""
        Compares two sets of basepairs, returns the counts of new and dropped basepairs in each category.
        """,
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("--baseline", type=str, required=False, help="Input CSV/Parquet file with the baseline data.")
    parser.add_argument("--target", type=str, required=False, help="Input CSV/Parquet file with the target data.")
    parser.add_argument("--comparison", type=str, required=False, help="Input CSV/Parquet file with the comparison already precomputed (has comparison_in_current and comparison_in_baseline columns).")
    parser.add_argument("--output", "-o", required=True, help="Output summary CSV file name.")
    parser.add_argument("--output-full-diff", "-O", required=False, help="Output file name of complete difference Parquet.")
    parser.add_argument("--boundaries", type=str, required=False, default=None, help="Return the results in the order of this boundaries file.")

    args = parser.parse_args()
    main(args)

