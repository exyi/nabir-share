import itertools
import subprocess
from typing import Any, Generator, Optional, Union
from matplotlib.axes import Axes
from matplotlib.figure import Figure
import matplotlib.patheffects
import polars as pl, numpy as np, numpy.typing as npt
import os, sys, math, re
import pairs
import pair_defs
from pair_defs import PairType
import pair_csv_parse
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
    # ("No filter ≤3 Å", (pl.col("resolution") <= 3)),
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
    def drop_columns(self, cols):
        ix = [ i for i, c in enumerate(self.columns) if c not in cols ]
        return self.select_columns(ix)
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
        ["hb_0_length", "hb_1_length", "hb_2_length", "hb_3_length", "hb_4_length"],
        # bin_width=0.05,
        pseudomin=2,
        pseudomax=6
    ),
    HistogramDef(
        "H-bond donor angle",
        "Angle (°)",
        ["hb_0_donor_angle", "hb_1_donor_angle", "hb_2_donor_angle", "hb_3_donor_angle", "hb_4_donor_angle"],
        # bin_width=2,
        pseudomin=0,
        pseudomax=360
    ),
    HistogramDef(
        "H-bond acceptor angle",
        "Angle (°)",
        ["hb_0_acceptor_angle", "hb_1_acceptor_angle", "hb_2_acceptor_angle", "hb_3_acceptor_angle", "hb_4_acceptor_angle"],
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

def is_symmetric_pair_type(pair_type: PairType):
    if not pair_type.swap_is_nop():
        return False

    hbonds = pair_defs.get_hbonds(pair_type)
    return all(
        pair_defs.hbond_swap_nucleotides(hb) in hbonds
        for hb in hbonds
    )

def get_label(col: str, pair_type: PairType, throw=True):
    hbonds = pair_defs.get_hbonds(pair_type, throw=throw)
    if not hbonds:
        assert not throw
        return None
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

def is_swapped(pair_type: PairType):
    symtype = pair_type.type in ["cWW", "tWW", "cHH", "tHH", "cSS", "tSS"]
    return pair_type.bases_str == "C-G" and symtype or pair_type.type.lower() not in pair_defs.pair_types
def format_pair_type(pair_type: PairType, is_dna = False, is_rna=False):
    pair_kind, pair_bases = pair_type.to_tuple()
    if is_dna == True:
        pair_bases = pair_bases.replace("U", "T")
    elif is_rna == True:
        pair_bases = pair_bases.replace("T", "U")
    if is_swapped(pair_type):
        return format_pair_type(pair_type.swap())
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

def get_bounds(dataframes: list[pl.DataFrame], pair_type: PairType, h: HistogramDef):
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
            return h.pseudomin or 0, h.pseudomax or 1
        all_datapoints = np.concatenate(datapoint_columns)
        mean = float(np.mean(all_datapoints))
        std = float(np.std(all_datapoints))
        pseudomin = max(mean - 3 * std, h.pseudomin or -math.inf)
        pseudomax = min(mean + 3 * std, h.pseudomax or math.inf)
        xmin = min(pseudomin, float(np.min([ np.quantile(c, 0.02) for c in datapoint_columns ])))
        xmax = max(pseudomax, float(np.max([ np.quantile(c, 0.98) for c in datapoint_columns ])))
        if xmin >= xmax - 0.0001:
            return min(xmin, h.pseudomin or 0), max(xmax, h.pseudomax or 1)
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

def get_hidden_columns(df: pl.DataFrame, pair_type: PairType):
    pt_hbonds = pair_defs.get_hbonds(pair_type, throw=True)
    hide_columns = [ i for i, hb in enumerate(pt_hbonds) if pair_defs.is_bond_hidden(pair_type, hb) ]
    # if len(hide_columns) == len(pt_hbonds):
    #     # all filtered out? fuckit
    #     visible_columns = h.columns

    columns_to_drop = df.limit(1).select(*(
        pl.col(f"^hb_{i}_.*$")
        for i in hide_columns
    )).columns
    print(f"{pair_type}: {hide_columns=} drop={columns_to_drop} hbonds={pt_hbonds}")
    return set(columns_to_drop)

def make_histogram_group(dataframes: list[pl.DataFrame], axes: list[Axes], titles: list[str], pair_type: PairType, h: HistogramDef):
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
                if i < len(legend) and j < len(legend) and legend[j] in renamed_columns.columns
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
        binses = np.arange(xmin, xmax, bin_width)
        sns.histplot(data=renamed_columns.to_pandas(),
                    #  binwidth=bin_width if len(renamed_columns) > 2 else None,
                    #  binwidth=bin_width,
                     bins=binses, # type:ignore
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
        #     if len(col) < 4:
        #         continue
        #     if len(col) > 30:
        #         col = col.filter((col > col.quantile(0.05)) & (col < col.quantile(0.95)))
        #     mean = float(col.mean())
        #     std = float(col.std())
        #     x = np.linspace(xmin, xmax, 100)
        #     y = scipy.stats.norm.pdf(x, mean, std) * len(renamed_columns) * bin_width
        #     ax.plot(x, y, color=f"C{i}", linestyle="--", linewidth=1)



def make_subplots(sp = subplots):
    fig, sp = plt.subplots(*sp)
    return fig, list(sp.reshape(-1))

def make_bond_pages(df: pl.DataFrame, outdir: str, pair_type: PairType, hs: list[HistogramDef], images = None, highlights: Optional[list[Optional[pl.DataFrame]]] = None, title_suffix = ""):
    hidden_bonds = get_hidden_columns(df, pair_type)
    if len(hidden_bonds) > 0:
        hs = [ h.drop_columns(hidden_bonds) for h in hs ]
    dataframes = [ df.filter(resolution_filter) for _, resolution_filter in resolutions ]
    # if sum(len(df) for df in dataframes) < 70:
    #     return
    pages: list[tuple[Figure, list[Axes]]] = [ make_subplots(subplots) for _ in resolutions ]
    titles = [ f"{format_pair_type(pair_type, is_dna=('DNA' in resolution_lbl))} {resolution_lbl}{title_suffix}" for (resolution_lbl, _), df in zip(resolutions, dataframes) ]
    print(titles)
    for p, title, df in zip(pages, titles, dataframes):
        fig, _ = p
        # fig.tight_layout(pad=3.0)
        fig.suptitle(title + f" ({len(df)}, class {determine_bp_class(df, pair_type)})")

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

            for ax, h in zip(p[1][1:], hs):
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

def make_resolution_comparison_page(df: pl.DataFrame, outdir: str, pair_type: PairType, h: HistogramDef, images = []):
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
    try:
        os.makedirs(outdir, exist_ok=True)
        pdf=os.path.join(outdir, title + ".pdf")
        fig.savefig(pdf, dpi=300)
        fig.savefig(os.path.join(outdir, title + ".png"))
        plt.close(fig)
        print(f"Wrote {pdf}")
        return pdf
    except Exception as e:
        print(f"Error writing {title}: {e}")
        raise e

def load_pair_table(file: str):
    df = pl.read_parquet(file) if file.endswith(".parquet") else pl.read_csv(file)
    df = pair_csv_parse.normalize_columns(df)

    if "coplanarity_angle" not in df.columns:
        # older file version
        df = df.with_columns(
            pl.col("^hb_\\d+_donor_angle$") / np.pi * 180,
            pl.col("^hb_\\d+_acceptor_angle$") / np.pi * 180,
        )
    return df

def infer_pair_type(filename: str):
    if m := re.match(r"^(n?[ct][HSW]{2}a?)-([AGCUT]-[AGCUT])\b", filename):
        return m.group(1), m.group(2)
    elif m := re.match(r"^([AGCUT]-[AGCUT])-(n?[ct][HSW]{2}a?)\b", filename):
        return m.group(2), m.group(1)
    else:
        return None
    
def tranpose_dict(d, columns):
    return {
        (c + "_" + k): v[i]
        for i, c in enumerate(columns)
        for k, v in d.items()
    }

def sample_for_kde(x: np.ndarray, threshold = 5_000):
    if len(x) <= threshold:
        return x
    else:
        return np.random.choice(x, threshold, replace=False)

def determine_bp_class(df: pl.DataFrame, pair_type: PairType, is_rna = None, throw=True):
    if is_rna is None:
        if len(df) == 0:
            is_rna = True
        else:
            is_rna = df['res1'].str.starts_with("D").not_().any() or df['res2'].str.starts_with("D").not_().any()
    assert isinstance(is_rna, bool)
    hbonds = pair_defs.get_hbonds(pair_type, throw=throw)
    if not hbonds:
        return None
    non_c_bonds = [ b for b in hbonds if not pair_defs.is_ch_bond(pair_type, b) ]
    good_base_bonds = [ b for b in hbonds if not pair_defs.is_bond_to_sugar(pair_type, b) and not pair_defs.is_ch_bond(pair_type, b) ]

    def unique_atoms1(bonds):
        return set(atom1 for _, atom1, _, _ in bonds)
    def unique_atoms2(bonds):# -> set[Any]:
        return set(atom2 for _, _, atom2, _ in bonds)
    
    def n_unique_atoms(bonds):
        return min(len(unique_atoms1(bonds)), len(unique_atoms2(bonds)))
    
    print(f"{pair_type} ~C={len(non_c_bonds)} good={len(good_base_bonds)} all={len(hbonds)} {non_c_bonds}")
    print(f"    {list(unique_atoms1(non_c_bonds))} {list(unique_atoms2(non_c_bonds))} | #uniq = {n_unique_atoms(good_base_bonds)} {n_unique_atoms(non_c_bonds)}")

    if n_unique_atoms(good_base_bonds) >= 2:
        return 1
    elif is_rna and n_unique_atoms(good_base_bonds) >= 1 and n_unique_atoms(non_c_bonds) >= 2:
        return 2
    else:
        return 3

def calculate_stats(df: pl.DataFrame, pair_type):
    if len(df) == 0:
        raise ValueError("No data")
    columns = df.select(pl.col("^hb_\\d+_(length|donor_angle|acceptor_angle)$")).columns
    # print(f"{columns=}")
    cdata = [ df[c].drop_nulls().to_numpy() for c in columns ]
    kdes = [ scipy.stats.gaussian_kde(sample_for_kde(c)) if len(c) > 5 else None for c in cdata ]
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
    # print(f"{nicest_basepairs=}")
    # print(f"{nicest_basepair=} LL={calc_datapoint_log_likelihood(df[nicest_basepair, :])} Σσ={score[nicest_basepair]} {next(df[nicest_basepair, columns].iter_rows())}")

    new_df_columns = {
        "mode_deviations": calc_datapoint_mode_deviations,
        "log_likelihood": calc_datapoint_log_likelihood,
    }

    result_stats = {
        "count": len(df),
        "bp_class": determine_bp_class(df, pair_type, throw=False),
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

def create_pair_image(row: pl.DataFrame, output_dir: str, pair_type: PairType) -> Optional[str]:
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

def reexport_df(df: pl.DataFrame, columns):
    df = df.with_columns(
        is_some_quality.alias("jirka_approves"),
        *[
            pl.Series(columns[col](df)).alias(col)
            for col in columns
        ]
    )
    df = df.drop([col for col in df.columns if re.match(r"[DR]NA-(0-1[.]8|1[.]8-3[.]5)(-r\d+)?", col)])
    # round float columns
    df = df.with_columns([
        pl.col(c).round_sig_figs(5).cast(pl.Float32).alias(c) for c in df.columns if df[c].dtype == pl.Float64 or df[c].dtype == pl.Float32
    ])
    return df

def enumerate_pair_types(files: list[str], include_nears: bool) -> Generator[tuple[PairType, pl.DataFrame], None, None]:
    assert len(files) > 0 and isinstance(files, list) and isinstance(files[0], str)
    for file in files:
        pair_type = infer_pair_type(os.path.basename(file))
        if pair_type is not None:
            yield PairType.from_tuple(pair_type), load_pair_table(file)
        else:
            df = load_pair_table(file)
            assert "type" in df.columns, f"{file} does not contain type column"
            df = df.with_columns(
                (pl.col("res1").replace(pairs.resname_map) + "-" +
                 pl.col("res2").replace(pairs.resname_map)
                ).alias("pair_bases")
            )
            groups = df.group_by("type", "pair_bases")
            print(f"{file}: {len(df)} rows, types: {list(sorted([ (k, len(gdf)) for k, gdf in groups ], key=lambda x: x[1], reverse=True))}")
            all_pairs_types = set(pair_defs.PairType.from_tuple(pt) for pt, _ in groups)
            for k, gdf in groups:
                k: Any
                pair_type = pair_defs.PairType.from_tuple(k)
                if len(set(pair_type.bases).difference(["A", "C", "G", "U", "T"])) > 0:
                    print(f"skipping weird bases: {pair_type}, count = {len(gdf)}")
                    continue
                if pair_type.is_swappable() and not pair_type.is_preferred_orientation() and pair_type.swap() in all_pairs_types:
                    print(f"skipping {pair_type} because it is redundant")
                    continue
                if pair_type.type[1].islower() and pair_type.type[2].isupper() and pair_type.type[1] == pair_type.type[2].lower():
                    continue
                if pair_type.n and not include_nears:
                    continue
                yield pair_type, gdf

def save_statistics(all_statistics, output_dir):
    df = pl.DataFrame(all_statistics, infer_schema_length=100_000)
    df.write_csv(os.path.join(output_dir, "statistics.csv"))
    bond_count = 10
    pt_family_dict = { pt: ix + 1 for ix, pt in enumerate(pair_defs.pair_types) }
    df2 = pl.concat([
        df.select(
            pl.col("bp_class").alias("Class"),
            pl.col("pair_type").str.to_lowercase().replace(pt_family_dict, default=-1).alias("Family"),
            pl.col("pair_type").alias("LW pair type"),
            pl.col("pair").alias("Pair bases"),
            pl.col("pair").str.split("-").map_elements(lambda x: x[0]).alias("Base 1"),
            pl.col("pair").str.split("-").map_elements(lambda x: x[1]).alias("Base 2"),
            pl.col("count").alias("Count"),
            pl.col("resolution_cutoff").alias("Resolution cutoff"),
            pl.lit(i).alias("hb_ix"),
            pl.col(f"hb_{i}_label").alias("H-bond Atoms"),
            pl.col(f"hb_{i}_length_mode").cast(pl.Float64).alias("Mode Distance"),
            pl.col(f"hb_{i}_length_median").cast(pl.Float64).alias("Median Distance"),
            pl.col(f"hb_{i}_length_mean").cast(pl.Float64).alias("Mean Distance"),
            pl.col(f"hb_{i}_length_std").cast(pl.Float64).alias("Std Distance"),
            pl.col(f"hb_{i}_donor_angle_mode").cast(pl.Float64).alias("Mode Donor Angle"),
            pl.col(f"hb_{i}_donor_angle_median").cast(pl.Float64).alias("Median Donor Angle"),
            pl.col(f"hb_{i}_donor_angle_mean").cast(pl.Float64).alias("Mean Donor Angle"),
            pl.col(f"hb_{i}_donor_angle_std").cast(pl.Float64).alias("Std Donor Angle"),
            pl.col(f"hb_{i}_acceptor_angle_mode").cast(pl.Float64).alias("Mode Acceptor Angle"),
            pl.col(f"hb_{i}_acceptor_angle_median").cast(pl.Float64).alias("Median Acceptor Angle"),
            pl.col(f"hb_{i}_acceptor_angle_mean").cast(pl.Float64).alias("Mean Acceptor Angle"),
            pl.col(f"hb_{i}_acceptor_angle_std").cast(pl.Float64).alias("Std Acceptor Angle"),
        )
        for i in range(bond_count)
        if f"hb_{i}_label" in df.columns
    ])
    df2 = df2.filter(pl.col("Mean Distance").is_not_null())
    hidden_col = []
    for pt, b, ix in zip(df2["LW pair type"], df2["Pair bases"], df2["hb_ix"]):
        hbonds = pair_defs.get_hbonds((pt, b))
        hidden = False
        if ix < len(hbonds):
            hidden = pair_defs.is_bond_hidden((pt, b), hbonds[ix])
        else:
            print(f"WARNING: {pt} {b} has only {len(hbonds)} bonds, but {ix} is requested ({hbonds})")

        hidden_col.append(hidden)

    df2 = df2.with_columns(
        pl.Series("hidden", hidden_col, dtype=pl.Boolean)
    )
    df2 = df2.sort([ "Class", "Family", "Base 1", "Base 2", "hb_ix" ])
    print("Wrote", os.path.join(output_dir, "statistics.csv"), "and", os.path.join(output_dir, "statistics2.csv"))
    df2.write_csv(os.path.join(output_dir, "statistics2.csv"))

    import xlsxwriter
    xlsx = os.path.join(output_dir, "statistics2.xlsx")
    with xlsxwriter.Workbook(xlsx) as workbook:

        df2.write_excel(workbook, worksheet="All", dtype_formats={ pl.Float64: "0.00" }, hidden_columns=["hb_ix", "Pair bases"])
        resolutions = df2["Resolution cutoff"].unique().to_list()
        for r in resolutions:
            df2.filter(pl.col("Resolution cutoff") == r)\
                .write_excel(workbook, worksheet=f"{r}", dtype_formats={ pl.Float64: "0.00" }, hidden_columns=["hb_ix", "Pair bases", "Resolution cutoff"])


def main(argv):
    import argparse
    parser = argparse.ArgumentParser(description="Generate histogram plots")
    parser.add_argument("input_file", help="Input file", nargs="+")
    parser.add_argument("--residue-directory", required=True)
    parser.add_argument("--reexport", default='none', choices=['none', 'partitioned'], help="Write out parquet files with calculated statistics columns (log likelihood, mode deviations)")
    parser.add_argument("--include-nears", default=False, action="store_true", help="If FR3D is run in basepair_detailed mode, it reports near basepairs (denoted as ncWW). By default, we ignore them, but this option includes them in the output.")
    parser.add_argument("--output-dir", "-o", required=True, help="Output directory")
    args = parser.parse_args(argv)

    residue_lists = residue_filter.read_id_lists(args.residue_directory)

    results = []

    all_statistics = []

    for pair_type, df in enumerate_pair_types(args.input_file, args.include_nears):
        df = residue_filter.add_res_filter_columns(df, residue_lists)
        print(f"{pair_type}: total count = {len(df)}, quality count = {len(df.filter(is_some_quality))}")
        print(df.select(pl.col("^hb_\\d+_length$"), pl.col("resolution"), is_some_quality.alias("some_quality")).describe())

        # good_bonds = [ i for i, bond in enumerate(pair_defs.get_hbonds(pair_type, throw=False) or []) if not pair_defs.is_bond_hidden(pair_type, bond) ]
        # df = df.filter(pl.all_horizontal(pl.lit(True), *[
        #     pl.col(f"hb_{i}_length") <= 3.1
        #     for i in good_bonds
        # ]))

        dff = None
        stat_columns = {
            "mode_deviations": lambda df: [ None ] * len(df),
            "log_likelihood": lambda df: [ None ] * len(df),
        }
        nicest_bps: Optional[list[int]] = None
        statistics = []
        for resolution_label, resolution_filter in {
            # "unfiltered": True,
            "1.8 Å": (pl.col("resolution") <= 1.8) & is_some_quality,
            "2.5 Å": (pl.col("resolution") <= 2.5) & is_some_quality,
            "RNA 3.0 Å": (pl.col("resolution") <= 3.0) & is_some_quality & is_rna,
            "DNA 3.0": (pl.col("resolution") <= 3.0) & is_some_quality & is_dna,
            "3.0 Å": (pl.col("resolution") <= 3.0) & is_some_quality,
        }.items():
            dff = df.filter(resolution_filter)\
                    .filter(pl.any_horizontal(pl.col("^hb_\\d+_length$").is_not_null()))
            if len(dff) == 0:
                continue
            stat_columns, stats = calculate_stats(dff, pair_type)
            statistics.append({
                "pair": pair_type.bases_str,
                "pair_type": pair_type.full_type,
                "resolution_cutoff": resolution_label,
                **stats,
            })
            if "nicest_bp_indices" in statistics[-1]:
                nicest_bps = statistics[-1]["nicest_bp_indices"]
                del statistics[-1]["nicest_bp_indices"]
            print(f"{pair_type} {resolution_label}: {len(dff)}/{len(df)} ")
        if dff is None or stat_columns is None:
            print(f"WARNING: No data in {pair_type} ({len(df)=}, len(filtered)={len(df.filter(is_some_quality))}")
            output_files = []
        elif not pair_defs.get_hbonds(pair_type, throw=False):
            print(f"WARNING: No hbonds for {pair_type}")
            output_files = []
        else:

            print("nicest_bps:", nicest_bps, "out of", len(dff) if dff is not None else 0)
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
        if args.reexport == "partitioned":
            reexport_df(df, stat_columns or []).write_parquet(os.path.join(args.output_dir, f"{pair_type}.parquet"))
            reexport_df(df.filter(is_some_quality), stat_columns or []).write_parquet(os.path.join(args.output_dir, f"{pair_type}-filtered.parquet"))
        results.append({
            # "input_file": file,
            "pair_type": pair_type.to_tuple(),
            "count": len(df),
            "high_quality": len(df.filter(is_high_quality)),
            "med_quality": len(df.filter(is_med_quality)),
            "score": len(df.filter(is_high_quality)) + len(df.filter(is_med_quality)) / 100,
            "files": output_files,
            "bp_class": statistics[-1]["bp_class"] if len(statistics) > 0 else determine_bp_class(df, pair_type, throw=False),
            "statistics": statistics,
            "labels": [
                get_label(f"hb_{i}_length", pair_type, throw=False) for i in range(3)
            ],
            "atoms": pair_defs.get_hbonds(pair_type, throw=False),
        })

    # results.sort(key=lambda r: r["score"], reverse=True)
    # results.sort(key=lambda r: r["pair_type"])
    results.sort(key=lambda r: (r["bp_class"] or 5, pair_defs.PairType.from_tuple(r["pair_type"])))
    output_files = [ f for r in results for f in r["files"] ]

    subprocess.run(["gs", "-dBATCH", "-dNOPAUSE", "-q", "-sDEVICE=pdfwrite", "-dPDFSETTINGS=/prepress", f"-sOutputFile={os.path.join(args.output_dir, 'hbonds-merged.pdf')}", *output_files])
    print("Wrote", os.path.join(args.output_dir, 'hbonds-merged.pdf'))
    save_statistics(all_statistics, args.output_dir)

    with open(os.path.join(args.output_dir, "output.json"), "w") as f:
        import json
        json.dump(results, f, indent=4)

if __name__ == "__main__":
    main(sys.argv[1:])

