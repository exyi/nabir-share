import itertools
import subprocess
from typing import Optional
import polars as pl, numpy as np
import os, sys, math, re
import pairs
import seaborn as sns
import matplotlib.pyplot as plt
import residue_filter
from dataclasses import dataclass


bins_per_width = 80
hist_kde = True

is_high_quality = pl.col("RNA-0-1.8") | pl.col("DNA-0-1.8")
is_some_quality = pl.col("RNA-0-1.8") | pl.col("DNA-0-1.8") | pl.col("DNA-1.8-3.5") | pl.col("RNA-1.8-3.5")
is_med_quality = pl.col("RNA-1.8-3.5") | pl.col("DNA-1.8-3.5")
is_dna = pl.col("res1").str.starts_with("D") | pl.col("res2").str.starts_with("D")
is_rna = pl.col("res1").str.starts_with("D").is_not() | pl.col("res2").str.starts_with("D").is_not()

plt.rcParams["figure.figsize"] = (16, 9)
subplots = (2, 2)
# subplots = None
# plt.rcParams["figure.figsize"] = (10, 6)

resolutions = [
    ("RNA ≤1.8 Å", is_high_quality & is_rna),
    ("RNA 1.8 Å - 3.5 Å", is_med_quality & is_rna),
    ("DNA ≤1.8 Å", is_high_quality & is_rna.is_not()),
    ("DNA 1.8 Å - 3.5 Å", is_med_quality & is_rna.is_not()),
    # ("≤1.8 Å", is_some_quality),
    # ("1.8 Å - 3.5 Å", is_med_quality),
]

@dataclass
class HistogramDef:
    title: str
    axis_label: str
    columns: list[str]
    legend: Optional[list[str]] = None
    bin_width: Optional[float] = None
    min: Optional[float] = None
    max: Optional[float] = None
    pseudomin: Optional[float] = None
    pseudomax: Optional[float] = None

histogram_defs = [
    HistogramDef(
        "H-bond length",
        "Distance (Å)",
        ["hb_0_length", "hb_1_length", "hb_2_length"],
        pseudomin=2,
        pseudomax=6
    ),
    HistogramDef(
        "H-bond donor angle",
        "Angle (°)",
        ["hb_0_donor_angle", "hb_1_donor_angle", "hb_2_donor_angle"],
        pseudomin=0,
        pseudomax=360
    ),
    HistogramDef(
        "H-bond acceptor angle",
        "Angle (°)",
        ["hb_0_acceptor_angle", "hb_1_acceptor_angle", "hb_2_acceptor_angle"],
        pseudomin=0,
        pseudomax=360
    )
]

def format_angle_label(atoms: tuple[str, str, str], swap=False):
    first = "B" if swap else "A"
    residue0 = atoms[0][0]
    if residue0 != first:
        atoms = tuple(reversed(atoms))

    return "".join([
        atoms[0][1:],
        "-" if atoms[1][0] == first else " · · · ",
        atoms[1][1:],
        " · · · " if atoms[1][0] == first else "-",
        atoms[2][1:]
    ])
def format_length_label(atoms: tuple[str, str], swap=False):
    a,b=atoms
    first = "B" if swap else "A"
    residue0 = atoms[0][0]
    if residue0 != first:
        a,b = b,a

    return f"{a[1:]} · · · {b[1:]}"

def get_label(col: str, pair_type: tuple[str, str]):
    hbonds = pairs.hbonding_atoms[pair_type]
    swap = is_swapped(pair_type[1])
    
    if col == "hb_0_donor_angle":
        return format_angle_label(hbonds[0][:3], swap=swap)
    elif col == "hb_1_donor_angle":
        return format_angle_label(hbonds[1][:3], swap=swap)
    elif col == "hb_2_donor_angle" and len(hbonds) > 2:
        return format_angle_label(hbonds[2][:3], swap=swap)
    elif col == "hb_0_acceptor_angle":
        return format_angle_label(hbonds[0][1:], swap=swap)
    elif col == "hb_1_acceptor_angle":
        return format_angle_label(hbonds[1][1:], swap=swap)
    elif col == "hb_2_acceptor_angle" and len(hbonds) > 2:
        return format_angle_label(hbonds[2][1:], swap=swap)
    elif col == "hb_0_length":
        return format_length_label(hbonds[0][1:3], swap=swap)
    elif col == "hb_1_length":
        return format_length_label(hbonds[1][1:3], swap=swap)
    elif col == "hb_2_length" and len(hbonds) > 2:
        return format_length_label(hbonds[2][1:3], swap=swap)
    
def is_swapped(pair_bases):
    return pair_bases == "C-G"
def format_pair_type(pair_type):
    pair_kind, pair_bases = pair_type
    if is_swapped(pair_bases):
        assert len(pair_kind) == 3
        return format_pair_type((f"{pair_kind[0]}{pair_kind[2]}{pair_kind[1]}", "-".join(reversed(pair_bases.split("-")))))
    elif len(pair_bases) == 3 and pair_bases[1] == '-':
        return pair_kind + " " + pair_bases.replace("-", "")
    else:
        return pair_kind + " " + pair_bases

def make_histogram(df: pl.DataFrame, outdir: str, pair_type: tuple[str, str], h: HistogramDef):
    title = f"{format_pair_type(pair_type)} {h.title}"
    if h.min is not None:
        assert h.max is not None
        xmin = h.min
        xmax = h.max
    else:
        datapoint_columns = [
            df.filter(f)[col].drop_nulls().to_numpy()
            for _, f in resolutions
            for col in h.columns
        ]
        datapoint_columns = [ c for c in datapoint_columns if len(c) > 0]
        if len(datapoint_columns) == 0:
            raise ValueError(f"No datapoints for {title}")
        all_datapoints = np.concatenate(datapoint_columns)
        mean = float(np.mean(all_datapoints))
        std = float(np.std(all_datapoints))
        pseudomin = max(mean - 3 * std, h.pseudomin or -math.inf)
        pseudomax = min(mean + 3 * std, h.pseudomax or math.inf)
        xmin = min(pseudomin, float(np.min([ np.quantile(c, 0.02) for c in datapoint_columns ])))
        xmax = max(pseudomax, float(np.max([ np.quantile(c, 0.98) for c in datapoint_columns ])))

    if h.bin_width is not None:
        bin_width = h.bin_width
    else:
        bin_width = (xmax - xmin) / bins_per_width

    if h.legend is not None:
        assert len(h.legend) == len(h.columns)
        legend = h.legend
    else:
        legend = [ get_label(col, pair_type) or "" for col in h.columns ]

    if subplots:
        main_fig, axes = plt.subplots(*subplots)
        main_fig.tight_layout(pad=3.0)
        axes = axes.reshape(len(resolutions))
    else:
        main_fig = None
        axes = [ None ] * len(resolutions)

    for ax_i, (resolution_lbl, resolution_filter) in enumerate(resolutions):
        dff = df.filter(resolution_filter)
        plot: plt.Axes = axes[ax_i] or plt.gca()
        plot.set(xlabel=h.axis_label,
                 title=f"{resolution_lbl}" if subplots else f"{title} {resolution_lbl}"
            )
        plot.set_xlim(xmin, xmax)
        nn_columns = [ c for c in h.columns if len(dff[c].drop_nulls()) > 0 ]

        if len(dff) < 1 or len(nn_columns) == 0:
            continue

        dffs = dff[nn_columns]
        print(bin_width, xmin, xmax, len(dffs), resolution_lbl, title)
        renamed_columns = dffs.select(*[
            pl.col(c).alias(l)
            for c, l in zip(h.columns, legend)
            if c in nn_columns
        ])
        sns.histplot(data=renamed_columns.to_pandas(), binwidth=bin_width, kde=(hist_kde and len(dffs) >= 5), legend=True, ax=plot)
        
        print([ l for l, c in zip(legend, h.columns) if c in nn_columns ])
        # plot.legend([ l for l, c in zip(legend, h.columns) if c in nn_columns ])
        if not subplots:
            yield save(plot.get_title(), outdir)
    
    if subplots:
        assert main_fig is not None
        main_fig.suptitle(title)
        yield save(title, outdir)


def save(title, outdir):
    os.makedirs(outdir, exist_ok=True)
    plt.savefig(os.path.join(outdir, title + ".pdf"))
    plt.savefig(os.path.join(outdir, title + ".png"))
    plt.close()
    return os.path.join(outdir, title + ".pdf")

def load_pair_table(file: str):
    df = pl.read_parquet(file) if file.endswith(".parquet") else pl.read_csv(file)

    df = df.with_columns(
        pl.col("hb_0_donor_angle") / np.pi * 180,
        pl.col("hb_1_donor_angle") / np.pi * 180,
        pl.col("hb_2_donor_angle") / np.pi * 180,
        pl.col("hb_0_acceptor_angle") / np.pi * 180,
        pl.col("hb_1_acceptor_angle") / np.pi * 180,
        pl.col("hb_2_acceptor_angle") / np.pi * 180,
    )
    return df

def infer_pair_type(filename: str):
    if m := re.match(r"^([ct][HSW]{2})-([AGCUT]-[AGCUT])\b", filename):
        return m.group(1), m.group(2)
    elif m := re.match(r"^([AGCUT]-[AGCUT])-([ct][HSW]{2})\b", filename):
        return m.group(2), m.group(1)
    else:
        raise ValueError(f"Unknown pair type: {filename}")

def main(argv):
    import argparse
    parser = argparse.ArgumentParser(description="Generate histogram plots")
    parser.add_argument("input_file", help="Input file", nargs="+")
    parser.add_argument("--residue-directory", required=True)
    parser.add_argument("--output-dir", "-o", required=True, help="Output directory")
    parser.add_argument("--pairing-type", default=None, type=str, help="For example tWS-A-G")
    args = parser.parse_args(argv)

    residue_lists = residue_filter.read_id_lists(args.residue_directory)

    results = []

    for file in args.input_file:
        pair_type = infer_pair_type(args.pairing_type or os.path.basename(file))
        df = load_pair_table(file)
        df = residue_filter.add_res_filter_columns(df, residue_lists)
        print(file)
        print(df[["hb_0_length", "hb_1_length", "hb_2_length"]].describe())

        output_files = [
            f
            for h in histogram_defs
            for f in make_histogram(df, args.output_dir, pair_type, h)
        ]
        results.append({
            "count": len(df),
            "score": len(df.filter(is_high_quality)) + len(df.filter(is_med_quality)) / 100,
            "files": output_files,
        })

    results.sort(key=lambda r: r["score"], reverse=True)
    output_files = [ f for r in results for f in r["files"] ]

    subprocess.run(["gs", "-dBATCH", "-dNOPAUSE", "-q", "-sDEVICE=pdfwrite", "-dPDFSETTINGS=/prepress", f"-sOutputFile={os.path.join(args.output_dir, 'hbonds-merged.pdf')}", *output_files])
    print("Wrote", os.path.join(args.output_dir, 'hbonds-merged.pdf'))

if __name__ == "__main__":
    main(sys.argv[1:])

