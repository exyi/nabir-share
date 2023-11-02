import itertools
import subprocess
from typing import Any, Optional, Union
from matplotlib.axes import Axes
from matplotlib.figure import Figure
import matplotlib.patheffects
import polars as pl, numpy as np, numpy.typing as npt
import os, sys, math, re
import pairs
import pair_defs
import seaborn as sns
import scipy.stats
import matplotlib.pyplot as plt
import residue_filter
from dataclasses import dataclass
import dataclasses

bins_per_width = 50
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
    ("DNA ≤3 Å", is_some_quality & is_rna.is_not() & (pl.col("resolution") <= 3)),
    ("RNA ≤3 Å", is_some_quality & is_rna & (pl.col("resolution") <= 3)),
    # ("DNA >3 Å", is_rna.is_not() & (pl.col("resolution") > 3)),
    # ("RNA >3 Å", is_rna & (pl.col("resolution") > 3)),
    # ("DNA ≤1.8 Å", is_high_quality & is_rna.is_not()),
    # ("DNA 1.8 Å - 3.5 Å", is_med_quality & is_rna.is_not()),
    # ("≤1.8 Å", is_high_quality),
    # ("1.8 Å - 3.5 Å", is_med_quality),
]

@dataclass(frozen=True)
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

    def copy(self, **kwargs):
        return dataclasses.replace(self, **kwargs)
    def select_columns(self, ix: Union[int, list[int]]):
        if isinstance(ix, int):
            ix = [ix]
        columns = [ self.columns[i] for i in ix ]
        legend = [ self.legend[i] for i in ix ] if self.legend is not None else None
        # title = f"{self.title} #{','.join(str(i+1) for i in ix)}"
        return self.copy(columns=columns, legend=legend)


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
        atoms = (atoms[2], atoms[1], atoms[0])

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

def is_symmetric_pair_type(pair_type: tuple[str, str]):
    if pair_type[1] != pair_type[1][::-1] or pair_type[0][1] != pair_type[0][2]:
        return False

    hbonds = pair_defs.get_hbonds(pair_type)
    return all(
        pair_defs.hbond_swap_nucleotides(hb) in hbonds
        for hb in hbonds
    )

def get_label(col: str, pair_type: tuple[str, str]):
    hbonds = pair_defs.get_hbonds(pair_type)
    swap = is_swapped(pair_type)
    
    if (m := re.match("^hb_(\\d+)_donor_angle", col)):
        ix = int(m.group(1))
        if len(hbonds) <= ix: return None
        return format_angle_label(hbonds[ix][:3], swap=swap)
    elif (m := re.match("^hb_(\\d+)_acceptor_angle", col)):
        ix = int(m.group(1))
        if len(hbonds) <= ix: return None
        return format_angle_label(hbonds[ix][1:], swap=swap)
    elif (m := re.match("^hb_(\\d+)_length", col)):
        ix = int(m.group(1))
        if len(hbonds) <= ix: return None
        return format_length_label(hbonds[ix][1:3], swap=swap)

def is_swapped(pair_type: tuple[str, str]):
    symtype = pair_type[0] in ["cWW", "tWW", "cHH", "tHH", "cSS", "tSS"]
    return pair_type[1] == "C-G" and symtype or pair_type[0] not in pair_defs.pair_types
def format_pair_type(pair_type, is_dna = False, is_rna=False):
    pair_kind, pair_bases = pair_type
    if is_dna == True:
        pair_bases = pair_bases.replace("U", "T")
    elif is_rna == True:
        pair_bases = pair_bases.replace("T", "U")
    if is_swapped(pair_type):
        assert len(pair_kind) == 3
        return format_pair_type((f"{pair_kind[0]}{pair_kind[2]}{pair_kind[1]}", "-".join(reversed(pair_bases.split("-")))))
    elif len(pair_bases) == 3 and pair_bases[1] == '-':
        return pair_kind + " " + pair_bases.replace("-", "")
    else:
        return pair_kind + " " + pair_bases
    
def crop_image(img: np.ndarray, padding = (0, 0, 0, 0)):
    if len(img.shape) == 2:
        img = img.reshape((*img.shape, 1))
    xbitmap = np.any(img != 0, axis=(1, 2))
    if not np.any(xbitmap):
        return img
    xlim = (np.argmax(xbitmap), len(xbitmap) - np.argmax(xbitmap[::-1]))
    ybitmap = np.any(img != 0, axis=(0, 2))
    ylim = (np.argmax(ybitmap), len(ybitmap) - np.argmax(ybitmap[::-1]))

    xlim = (max(0, xlim[0] - padding[0]), min(img.shape[0], xlim[1] + padding[1]))
    ylim = (max(0, ylim[0] - padding[2]), min(img.shape[1], ylim[1] + padding[3]))

    return img[xlim[0]:xlim[1], ylim[0]:ylim[1], :]

def get_bounds(dataframes: list[pl.DataFrame], pair_type: tuple[str, str], h: HistogramDef):
    if h.min is not None:
        assert h.max is not None
        return h.min, h.max
    else:
        datapoint_columns = [
            df[col].drop_nulls().to_numpy()
            for df in dataframes
            for col in h.columns if col in df.columns
        ]
        datapoint_columns = [ c for c in datapoint_columns if len(c) > 0]
        if len(datapoint_columns) == 0:
            # raise ValueError(f"No datapoints for histogram {h.title} {pair_type}")
            return 0, 1
        all_datapoints = np.concatenate(datapoint_columns)
        mean = float(np.mean(all_datapoints))
        std = float(np.std(all_datapoints))
        pseudomin = max(mean - 3 * std, h.pseudomin or -math.inf)
        pseudomax = min(mean + 3 * std, h.pseudomax or math.inf)
        xmin = min(pseudomin, float(np.min([ np.quantile(c, 0.02) for c in datapoint_columns ])))
        xmax = max(pseudomax, float(np.max([ np.quantile(c, 0.98) for c in datapoint_columns ])))
        return xmin, xmax

def get_histogram_ticksize(max, max_ticks = 8):
    def it(max):
        yield 1
        yield 2
        yield 5
        yield 10
        for x in it(max/10):
            yield x*10
    for ticksize in it(max):
        if ticksize * max_ticks > max:
            return ticksize

def make_histogram_group(dataframes: list[pl.DataFrame], axes: list[plt.Axes], titles: list[str], pair_type: tuple[str, str], h: HistogramDef):
    xmin, xmax = get_bounds(dataframes, pair_type, h)

    if h.bin_width is not None:
        bin_width = h.bin_width
    else:
        bin_width = (xmax - xmin) / bins_per_width

    if h.legend is not None:
        assert len(h.legend) == len(h.columns)
        legend = h.legend
    else:
        legend = [ get_label(col, pair_type) or "" for col in h.columns ]
        legend = [ l for l in legend if l ]

    is_symmetric = is_symmetric_pair_type(pair_type)

    for df, ax, title in zip(dataframes, axes, titles):
        ax.set(xlabel=h.axis_label, title=title)
        ax.set_xlim(xmin, xmax)
        nn_columns = [ c for c in h.columns if c in df.columns and len(df[c].drop_nulls()) > 0 ]

        if len(df) < 1 or len(nn_columns) == 0:
            continue

        dfs = df[nn_columns]
        renamed_columns = dfs.select(*[
            pl.col(c).alias(l)
            for c, l in zip(h.columns, legend)
            if c in nn_columns
        ])
        if is_symmetric:
            print(f"{pair_type} is symetric")
            # merge symetric bonds
            all_hbons = pair_defs.get_hbonds(pair_type)
            symetric_bonds = list(set(
                tuple(sorted((i, all_hbons.index(pair_defs.hbond_swap_nucleotides(hb)))))
                for i, hb in enumerate(all_hbons)
            ))
            assert (0, 0) not in symetric_bonds and (1, 1) not in symetric_bonds and (2, 2) not in symetric_bonds
            renamed_columns_ = renamed_columns.select(
                pl.col(legend[j]).alias(legend[i])
                for i, j in symetric_bonds
                if i < len(legend) and j < len(legend)
            )
            if len(renamed_columns_.columns) == len(renamed_columns.columns):
                print("Merging symetric bonds: ", symetric_bonds)
                renamed_columns = renamed_columns.select(pl.col(legend[i]) for i, _ in symetric_bonds)
                renamed_columns = pl.concat([ renamed_columns, renamed_columns_ ])

        if len(renamed_columns.columns) == 0 or len(renamed_columns) == 0:
            print(f"WARNING: no columns left after merging ({title})")
            continue
        if len(renamed_columns) == 1:
            print(renamed_columns)

        print(bin_width, xmin, xmax, len(renamed_columns), len(renamed_columns.columns), renamed_columns.null_count().to_dicts()[0], title)
        # binses = np.arange(xmin, xmax, bin_width)
        sns.histplot(data=renamed_columns.to_pandas(),
                    #  binwidth=bin_width if len(renamed_columns) > 2 else None,
                     kde=(hist_kde and len(dfs) >= 5),
                     legend=True,
                     ax=ax)
        ymax = ax.get_ylim()[1]
        ax.set_yticks(np.arange(0, ymax, step=get_histogram_ticksize(ymax)))
        if hist_kde:
            for line in ax.lines:
                line: Any
                x = line.get_xdata()
                y = line.get_ydata()
                peak = np.argmax(y)
                # circle marker for the peak
                ax.plot([ x[peak] ], [ y[peak] ], marker="o", color=line.get_color())
                # text label for the peak
                peak_fmt = f"{x[peak]:.0f}°" if "angle" in title.lower() else f"{x[peak]:.2f}"
                ax.annotate(peak_fmt, (x[peak], y[peak]), xytext=(0, 5), textcoords="offset points", ha='center', va='bottom', color=line.get_color(), fontsize=8, path_effects=[
                    matplotlib.patheffects.withStroke(linewidth=3, foreground="white") # text outline
                ])

                curve_steps = np.append(x[1:] - x[:-1], [0])
                curve_area_total = np.sum(y * curve_steps)
                curve_area_cumsum = np.cumsum(y * curve_steps)
                quantiles = [ np.searchsorted(curve_area_cumsum, q * curve_area_total) for q in [ 0.05, 0.95 ] ]

                for q in quantiles:
                    ax.plot([ x[q] ], [ y[q] ], marker="|", markersize=10, color=line.get_color())

        
        # for i, colname in enumerate(renamed_columns.columns):

        #     # add normal distribution for comparison
        #     col = renamed_columns[colname]
        #     if len(col) > 10:
        #         col = col.filter((col > col.quantile(0.05)) & (col < col.quantile(0.95)))
        #     mean = float(col.mean())
        #     std = float(col.std())
        #     x = np.linspace(xmin, xmax, 100)
        #     y = scipy.stats.norm.pdf(x, mean, std) * len(renamed_columns) * bin_width
        #     ax.plot(x, y, color=f"C{i}", linestyle="--", linewidth=1)



def make_subplots(sp = subplots):
    fig, sp = plt.subplots(*sp)
    return fig, list(sp.reshape(-1))

def make_bond_pages(df: pl.DataFrame, outdir: str, pair_type: tuple[str, str], hs: list[HistogramDef], images = None, highlights: Optional[list[Optional[pl.DataFrame]]] = None, title_suffix = ""):

    dataframes = [ df.filter(resolution_filter) for _, resolution_filter in resolutions ]
    pages: list[tuple[Figure, list[Axes]]] = [ make_subplots(subplots) for _ in resolutions ]
    titles = [ f"{format_pair_type(pair_type, is_dna=('DNA' in resolution_lbl))} {resolution_lbl}{title_suffix}" for (resolution_lbl, _), df in zip(resolutions, dataframes) ]
    print(titles)
    for p, title, df in zip(pages, titles, dataframes):
        fig, _ = p
        # fig.tight_layout(pad=3.0)
        fig.suptitle(title + f" ({len(df)})")

    for i, h in enumerate(hs):
        make_histogram_group(dataframes, [ p[1][i+1] for p in pages ], [h.title] * len(dataframes), pair_type, h)
    if images is not None:
        for p, img, highlight in zip(pages, images, highlights or itertools.repeat(None)):
            if img is None:
                continue
            ax = p[1][0]
            img_data=crop_image(plt.imread(img), padding=(0, 30, 0, 0))
            print(f"image {img} {img_data.shape}")
            ax.imshow(img_data)
            # ax.annotate("bazmek", (0, 1))
            if highlight is not None:
                def fmt_ins(ins):
                    if not ins or ins == '\0' or ins == ' ':
                        return ""
                    else:
                        return ".ins" + ins
                def fmt_alt(alt):
                    if not alt or alt == '?':
                        return ""
                    else:
                        return ".alt" + alt
                chain1 = f"chain {highlight[0, 'chain1']} "
                chain2 = f"chain {highlight[0, 'chain2']} " if highlight[0, 'chain2'] != highlight[0, 'chain1'] else ""
                address = f"{highlight[0, 'pdbid']}: {chain1}{highlight[0, 'res1']}{highlight[0, 'nr1']}{fmt_ins(highlight[0, 'ins1'])}{fmt_alt(highlight[0, 'alt1'])} · · · {chain2}{highlight[0, 'res2']}{highlight[0, 'nr2']}{fmt_ins(highlight[0, 'ins2'])}{fmt_alt(highlight[0, 'alt2'])}"
                ax.text(0.5, 0, f"{address} ({highlight[0, 'resolution']:.1f} Å)", transform=ax.transAxes, horizontalalignment="center")
            # ax.legend("bazmek")
            ax.axis("off")
    if highlights is not None:
        for p, highlight in zip(pages, highlights):
            if highlight is None or len(highlight) == 0:
                continue

            for ax, h in zip(p[1][1:], histogram_defs):
                for col_i, col in enumerate(h.columns):
                    if col in highlight.columns and highlight[0, col] is not None:
                        ax.plot([ float(highlight[0, col]) ], [ 0 ], marker="o", color=f"C{col_i}")

    for p, title, df in zip(pages, titles, dataframes):
        if len(df) == 0:
            # make "NO DATA" page
            fig, ax = plt.subplots(1)
            fig.suptitle(title)
            ax.axis("off")
            ax.text(0.5, 0.5,'NO DATA',fontsize=30,horizontalalignment='center',verticalalignment='center',transform = ax.transAxes)
            yield save(fig, title, outdir)
        else:
            fig, axes = p
            # fig.tight_layout(pad=3.0)
            yield save(fig, title, outdir)

def make_resolution_comparison_page(df: pl.DataFrame, outdir: str, pair_type: tuple[str, str], h: HistogramDef, images = []):
    title = f"{format_pair_type(pair_type)} {h.title}"

    dataframes = [ df.filter(resolution_filter) for _, resolution_filter in resolutions ]

    if subplots:
        titles = [ f"{resolution_lbl}" for resolution_lbl, _ in resolutions ]
        main_fig, axes = plt.subplots(*subplots)
        main_fig.tight_layout(pad=3.0)
        axes = list(axes.reshape(-1))
        assert len(axes) == len(resolutions) + len(images)
        for ax_i, img in enumerate(images):
            ax_i += len(resolutions)
            axes[ax_i].imshow(plt.imread(img))
    else:
        titles = [ f"{title} {resolution_lbl}" for resolution_lbl, _ in resolutions ]
        main_fig = None
        axes = [ plt.gca() for _ in resolutions ]
    make_histogram_group(dataframes, axes, titles, pair_type, h)

    if subplots:
        assert main_fig is not None
        main_fig.suptitle(title)
        yield save(main_fig, title, outdir)
    else:
        for ax_i, (resolution_lbl, resolution_filter) in enumerate(resolutions):
            yield save(axes[ax_i].figure, axes[ax_i].get_title(), outdir) #type:ignore


def save(fig: Figure, title, outdir):
    os.makedirs(outdir, exist_ok=True)
    pdf=os.path.join(outdir, title + ".pdf")
    fig.savefig(pdf, dpi=300)
    fig.savefig(os.path.join(outdir, title + ".png"))
    plt.close(fig)
    print(f"Wrote {pdf}")
    return pdf

def load_pair_table(file: str):
    df = pl.read_parquet(file) if file.endswith(".parquet") else pl.read_csv(file)

    df = df.with_columns(
        pl.col("^hb_\\d+_donor_angle$") / np.pi * 180,
        pl.col("^hb_\\d+_acceptor_angle$") / np.pi * 180,
    )
    return df

def infer_pair_type(filename: str):
    if m := re.match(r"^([ct][HSW]{2})-([AGCUT]-[AGCUT])\b", filename):
        return m.group(1), m.group(2)
    elif m := re.match(r"^([AGCUT]-[AGCUT])-([ct][HSW]{2})\b", filename):
        return m.group(2), m.group(1)
    else:
        raise ValueError(f"Unknown pair type: {filename}")
    
def tranpose_dict(d, columns):
    return {
        (c + "_" + k): v[i]
        for i, c in enumerate(columns)
        for k, v in d.items()
    }

def calculate_stats(df: pl.DataFrame, pair_type):
    if len(df) == 0:
        raise ValueError("No data")
    columns = df.select(pl.col("^hb_\\d+_(length|donor_angle|acceptor_angle)$")).columns
    print(f"{columns=}")
    cdata = [ df[c].drop_nulls().to_numpy() for c in columns ]
    kdes = [ scipy.stats.gaussian_kde(c) if len(c) > 5 else None for c in cdata ]
    kde_modes = [ None if kde is None else c[np.argmax(kde.pdf(c))] for kde, c in zip(kdes, cdata) ]
    medians = [ np.median(c) if len(c) > 0 else None for c in cdata ]
    means = [ np.mean(c) if len(c) > 0 else None for c in cdata ]
    stds = [ np.std(c) if len(c) > 1 else None for c in cdata ]
    kde_mode_stds = [ None if kde_mode is None else np.sqrt(np.mean((c - kde_mode) ** 2)) for kde_mode, c in zip(kde_modes, cdata) ]
    def calc_datapoint_mode_deviations(df):
        return np.sum([
            ((df[c] - kde_mode) ** 2 / kde_mode_std ** 2).fill_null(10).to_numpy()
            for kde_mode, kde_mode_std, c in zip(kde_modes, kde_mode_stds, columns)
            if kde_mode is not None and kde_mode_std
        ] + [ [0] * len(df) ], axis=0)
    def calc_datapoint_log_likelihood(df):
        return np.sum([ kde.logpdf(df[c].fill_null(0).to_numpy()) if kde is not None else np.zeros(len(df)) for kde, c in zip(kdes, columns) ], axis=0)
    # print("datapoint_mode_deviations:", datapoint_mode_deviations)
    # print("datapoint_log_likelihood:", datapoint_log_likelihood)
    # nicest_basepair = int(np.argmax(datapoint_log_likelihood))
    score = -calc_datapoint_mode_deviations(df)
    min_score = np.min(score)
    score = score - min_score + 1
    nicest_basepair = int(np.argmax(score))
    nicest_basepairs = [
        int(np.argmax(np.concatenate([ [0.1], score * df.select(r.alias("x"))["x"].fill_null(False).to_numpy()]))) - 1
        for _, r in resolutions
    ]
    print(f"{nicest_basepairs=}")
    print(f"{nicest_basepair=} LL={calc_datapoint_log_likelihood(df[nicest_basepair, :])} Σσ={score[nicest_basepair]} {next(df[nicest_basepair, columns].iter_rows())}")

    new_df_columns = {
        "mode_deviations": calc_datapoint_mode_deviations,
        "log_likelihood": calc_datapoint_log_likelihood,
    }

    result_stats = {
        "count": len(df),
        "nicest_bp": str(next(df[nicest_basepair, ["pdbid", "model", "chain1", "res1", "nr1", "ins1", "alt1", "chain2", "res2", "nr2", "ins2", "alt2"]].iter_rows())),
        "nicest_bp_index": nicest_basepair,
        "nicest_bp_indices": nicest_basepairs,
        **{
            f"hb_{i}_label": get_label(f"hb_{i}_length", pair_type) for i in range(len(pair_defs.get_hbonds(pair_type)))
        },
        **tranpose_dict({
            "mode": kde_modes,
            "median": medians,
            "mean": means,
            "std": stds,
        }, columns)
    }

    return new_df_columns, result_stats

def create_pair_image(row: pl.DataFrame, output_dir: str, pair_type: tuple[str,str]) -> Optional[str]:
    if len(row) == 0:
        return None
    os.makedirs(os.path.join(output_dir, "img"), exist_ok=True)
    row.write_parquet(os.path.join(output_dir, "img", f"nicest.parquet"))
    pdbid = row["pdbid"]
    label_atoms = list(itertools.chain(*[ (x, y) for (_, x, y, _) in pair_defs.get_hbonds(pair_type) ]))
    command = [
        "pymol", "-cq",
        os.path.join(os.path.dirname(__file__), "gen_contact_images.py"),
        "--",
        os.path.join(output_dir, "img", f"nicest.parquet"),
        "--label-atoms", *label_atoms,
        f"--output-dir={os.path.join(output_dir, 'img')}",
    ]
    print(*command)
    p = subprocess.run(command, capture_output=True)
    if p.returncode != 0:
        print(p.stdout.decode('utf-8'))
        print(p.stderr.decode('utf-8'))
        raise ValueError(f"PyMOL failed with code {p.returncode}")
    output_lines = p.stdout.decode('utf-8').splitlines()
    for l in output_lines:
        if m := re.match(r"^Saved basepair image (.*)$", l):
            return m.group(1)
    print(*command)
    print(p.stdout.decode('utf-8'))
    print(p.stderr.decode('utf-8'))
    raise ValueError(f"Could not find PyMOL generated image file")

def reexport_df(df, columns):
    df = df.with_columns(
        is_some_quality.alias("jirka_approves"),
        *[
            pl.Series(columns[col](df)).alias(col)
            for col in columns
        ]
    )
    return df

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

    all_statistics = []

    for file in args.input_file:
        pair_type = infer_pair_type(args.pairing_type or os.path.basename(file))
        df = load_pair_table(file)
        df = residue_filter.add_res_filter_columns(df, residue_lists)
        print(file)
        print(df.select(pl.col("^hb_\\d+_length$")).describe())

        dff, stat_columns = None, None
        nicest_bps: Optional[list[int]] = None
        statistics = []
        for resolution_cutoff in [1.8, 2.2, 3.0]:
            dff = df.filter(is_some_quality)\
                    .filter(pl.col("resolution") <= resolution_cutoff)\
                    .filter(pl.any_horizontal(pl.col("^hb_\\d+_length$").is_not_null()))
            if len(dff) == 0:
                continue
            stat_columns, stats = calculate_stats(dff, pair_type)
            statistics.append({
                "pair": pair_type[1],
                "pair_type": pair_type[0],
                "resolution_cutoff": resolution_cutoff,
                **stats,
            })
            if "nicest_bp_indices" in statistics[-1]:
                nicest_bps = statistics[-1]["nicest_bp_indices"]
                del statistics[-1]["nicest_bp_indices"]
        if dff is None or stat_columns is None:
            print(f"WARNING: No data in {file}")
            continue

        print(nicest_bps, len(dff) if dff is not None else 0)
        # output_files = [
        #     f
        #     for h in histogram_defs
        #     for f in make_resolution_comparison_page(df, args.output_dir, pair_type, h, images= [ create_pair_image(df[nicest_bp], args.output_dir, pair_type) ] if nicest_bp is not None else [])
        # ]
        dna_rna_images = [ create_pair_image(dff[bp], args.output_dir, pair_type) if bp >= 0 else None for bp in nicest_bps ] * len(resolutions) if nicest_bps is not None else []
        dna_rna_highlights = [ dff[bp] if bp >= 0 else None for bp in nicest_bps ] if nicest_bps is not None else []
        output_files = [
            f for f in make_bond_pages(df, args.output_dir, pair_type, histogram_defs, images=dna_rna_images, highlights=dna_rna_highlights
            )
        ]
        # output_files = [
        #     f
        #     for column in [0, 1, 2]
        #     for f in make_bond_pages(df, args.output_dir, pair_type, [ h.select_columns(column) for h in histogram_defs], images=dna_rna_images, highlights=dna_rna_highlights, title_suffix=f" #{column}")
        # ]
        all_statistics.extend(statistics)
        reexport_df(df, stat_columns).write_parquet(os.path.join(args.output_dir, f"{pair_type[1]}-{pair_type[0]}.parquet"))
        results.append({
            "input_file": file,
            "pair_type": pair_type,
            "count": len(df),
            "score": len(df.filter(is_high_quality)) + len(df.filter(is_med_quality)) / 100,
            "files": output_files,
            "statistics": statistics,
            "labels": [
                get_label(f"hb_{i}_length", pair_type) for i in range(3)
            ],
            "atoms": pair_defs.get_hbonds(pair_type),
        })

    # results.sort(key=lambda r: r["score"], reverse=True)
    # results.sort(key=lambda r: r["pair_type"])
    results.sort(key=lambda r: (pair_defs.pair_types.index(r["pair_type"][0]), r["pair_type"][1][::-1]))
    output_files = [ f for r in results for f in r["files"] ]

    subprocess.run(["gs", "-dBATCH", "-dNOPAUSE", "-q", "-sDEVICE=pdfwrite", "-dPDFSETTINGS=/prepress", f"-sOutputFile={os.path.join(args.output_dir, 'hbonds-merged.pdf')}", *output_files])
    print("Wrote", os.path.join(args.output_dir, 'hbonds-merged.pdf'))
    statistics_df = pl.DataFrame(all_statistics)
    statistics_df.write_csv(os.path.join(args.output_dir, "statistics.csv"))
    print("Wrote", os.path.join(args.output_dir, "statistics.csv"))

    with open(os.path.join(args.output_dir, "output.json"), "w") as f:
        import json
        json.dump(results, f, indent=4)

if __name__ == "__main__":
    main(sys.argv[1:])

