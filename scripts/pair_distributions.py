#!/usr/bin/env python3
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
is_rna = pl.col("res1").str.starts_with("D").not_() | pl.col("res2").str.starts_with("D").not_()

plt.rcParams["figure.figsize"] = (16, 9)
subplots = (2, 2)
# subplots = None
# plt.rcParams["figure.figsize"] = (10, 6)

resolutions = [
    # ("DNA ≤3 Å", is_some_quality & is_rna.not_() & (pl.col("resolution") <= 3)),
    ("RNA ≤3 Å", is_some_quality & is_rna & (pl.col("resolution") <= 3)),
    # ("RNA ≤3 Å", is_some_quality & (pl.col("resolution") <= 3)),
    # ("No filter ≤3 Å", (pl.col("resolution") <= 3)),
    # ("DNA >3 Å", is_rna.not_() & (pl.col("resolution") > 3)),
    # ("RNA >3 Å", is_rna & (pl.col("resolution") > 3)),
    # ("DNA ≤1.8 Å", is_high_quality & is_rna.not_()),
    # ("DNA 1.8 Å - 3.5 Å", is_med_quality & is_rna.not_()),
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


hbond_histogram_defs = [
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

coplanarity_histogram_defs = [
    HistogramDef(
        "Yaw, pitch, left-to-right, N1-C1' reference frames",
        "Angle (°)",
        ["C1_C1_yaw1", "C1_C1_pitch1", "C1_C1_roll1"],
        legend=["Yaw 1", "Pitch 1", "Roll 1"],
        # bin_width=2,
        min=-180,
        max=180
    ),
    HistogramDef(
        "Euler angles, N1-C1' reference frames",
        "Angle (°)",
        ["C1_C1_euler_phi", "C1_C1_euler_theta", "C1_C1_euler_psi", "C1_C1_euler_phicospsi"],
        legend=["φ", "θ", "ψ", "φ-cos(θ)ψ"],
        # bin_width=0.05,
        min=-180,
        max=180
    ),
    HistogramDef(
        "Yaw, pitch, right-to-left, N1-C1' reference frames",
        "Angle (°)",
        ["C1_C1_yaw2", "C1_C1_pitch2", "C1_C1_roll2"],
        legend=["Yaw 2", "Pitch 2", "Roll 2"],
        min=-180,
        max=180
    ),
    HistogramDef(
        "RMSD of pairing edges",
        "Å",
        ["rmsd_edge1", "rmsd_edge2"],
        legend=["fit on left base", "fit on right base"],
    ),

]

coplanarity_histogram_defs2 = [
    HistogramDef(
        "H-bond angle to left plane",
        "Angle (°)",
        ["hb_0_OOPA1", "hb_1_OOPA1", "hb_2_OOPA1"],
        pseudomin=-180, pseudomax=180,
        min=-60, max=60
    ),
    HistogramDef(
        "H-bond angle to right plane",
        "Angle (°)",
        ["hb_0_OOPA2", "hb_1_OOPA2", "hb_2_OOPA2"],
        pseudomin=-180, pseudomax=180,
        min=-60, max=60
    ),
    # HistogramDef(
    #     "H-bond angle to donor plane",
    #     "Angle (°)",
    #     ["hb_0_donor_OOPA", "hb_1_donor_OOPA", "hb_2_donor_OOPA"],
    #     pseudomin=-180, pseudomax=180,
    #     min=-60, max=60
    # ),
    # HistogramDef(
    #     "H-bond angle to donor plane",
    #     "Angle (°)",
    #     ["hb_0_acceptor_OOPA", "hb_1_acceptor_OOPA", "hb_2_acceptor_OOPA"],
    #     pseudomin=-180, pseudomax=180,
    #     min=-60, max=60
    # ),
    HistogramDef(
        "Edge-to-plane distance",
        "Distance (Å)",
        ["coplanarity_shift1", "coplanarity_shift2"],
        legend=["Edge 2 - plane 1", "Edge 1 - Plane 2"],
        pseudomin=-6, pseudomax=6,
        min=-1.1, max=1.1
    ),
    HistogramDef(
        "Edge-to-plane angle",
        "Angle (°)",
        ["coplanarity_edge_angle1", "coplanarity_edge_angle2"],
        legend=["Edge 2 - plane 1", "Edge 1 - Plane 2"],
        pseudomin=-180, pseudomax=180,
        min=-60, max=60
    ),
]

coplanarity_histogram_defs_selection = [
    HistogramDef(
        "Yaw, pitch, left-to-right, N1-C1' reference frames",
        "Angle (°)",
        ["C1_C1_yaw1", "C1_C1_pitch1", "C1_C1_roll1"],
        legend=["Yaw 1", "Pitch 1", "Roll 1"],
        # bin_width=2,
        min=-180,
        max=180
    ),
    # HistogramDef(
    #     "H-bond angle to acceptor plane",
    #     "Angle (°)",
    #     ["hb_0_acceptor_OOPA", "hb_1_acceptor_OOPA", "hb_2_acceptor_OOPA"],
    #     pseudomin=-180, pseudomax=180,
    #     min=-60, max=60
    # ),

    HistogramDef(
        "Edge-to-plane angle",
        "Angle (°)",
        ["coplanarity_edge_angle1", "coplanarity_edge_angle2"],
        legend=["Edge 2 - plane 1", "Edge 1 - Plane 2"],
        pseudomin=-180, pseudomax=180,
        min=-60, max=60
    ),
    HistogramDef(
        "Yaw, pitch, right-to-left, N1-C1' reference frames",
        "Angle (°)",
        ["C1_C1_yaw2", "C1_C1_pitch2", "C1_C1_roll2"],
        legend=["Yaw 2", "Pitch 2", "Roll 2"],
        min=-180,
        max=180
    ),
    # HistogramDef(
    #     "H-bond angle to donor plane",
    #     "Angle (°)",
    #     ["hb_0_donor_OOPA", "hb_1_donor_OOPA", "hb_2_donor_OOPA"],
    #     pseudomin=-180, pseudomax=180,
    #     min=-60, max=60
    # ),

    # HistogramDef(
    #     "RMSD of pairing edges",
    #     "Å",
    #     ["rmsd_edge1", "rmsd_edge2"],
    #     legend=["fit on left base", "fit on right base"],
    # ),
]

rmsd_histogram_defs = [
    HistogramDef(
        "RMSD of C1'-N bonds",
        "Å",
        ["rmsd_C1N_frames1", "rmsd_C1N_frames2"],
        legend=["fit on left nucleotide", "fit on right nucleotide"],
    ),
    HistogramDef(
        "RMSD of pairing edges",
        "Å",
        ["rmsd_edge1", "rmsd_edge2"],
        legend=["fit on left base", "fit on right base"],
    ),
    HistogramDef(
        "RMSD of C1'-N bonds, ",
        "Å",
        ["rmsd_edge_C1N_frame"],
        legend=["fit on pairing edges"]
    ),
    HistogramDef(
        "Oveall basepair RMSD",
        "Å",
        ["rmsd_all_base"],
        legend=["fit on all base atoms"]
    ),
]

def is_angular_modular(column: str):
    return re.match(r".*(\b|_)(yaw|pitch|roll|euler_theta|euler_psi|euler_phi)\d*$", column) is not None

def angular_modulo(x):
    return (x + (180 + 2*360)) % 360 - 180
    #                  ^^^^ polars broken modulo workaround

def angular_modular_mean(data: Union[pl.Series, np.ndarray], throw_if_empty = False) -> Optional[float]:
    """
    >>> np.angle(np.mean(np.exp(1j * np.radians([-170, 175]))), deg=True)
    -177.5
    """
    if isinstance(data, pl.Series):
        data = data.drop_nulls().to_numpy()

    if len(data) == 0:
        if throw_if_empty:
            raise ValueError("Empty mean")
        else:
            return None
    
    return float(np.angle(np.mean(np.exp(1j * np.radians(data))), deg=True))

def angular_modular_std(data: Union[pl.Series, np.ndarray], throw_if_empty = False, mean = None) -> Optional[float]:
    """
    """
    if isinstance(data, pl.Series):
        data = data.drop_nulls().to_numpy()

    if len(data) <= 1:
        if throw_if_empty:
            raise ValueError("STD on empty or singleton set")
        else:
            return None
    on_circle = np.exp(1j * np.radians(data))
    if mean is None:
        mean_complex = np.mean(on_circle)
    else:
        mean_complex = np.exp(1j * np.radians(mean))
    around_zero = np.angle(on_circle / mean_complex, deg=True) # multiplication is rotation, division is rotation back
    return math.sqrt(float(np.mean(around_zero ** 2)))

def angular_modular_minmax(data: Union[pl.Series, np.ndarray], throw_if_empty = False, percentiles = (0, 100)) -> Optional[tuple[float, float]]:
    """
    """
    if isinstance(data, pl.Series):
        data = data.drop_nulls().to_numpy()
    if len(data) <= 0:
        if throw_if_empty:
            raise ValueError("Empty set has no min/max")
        else:
            return None

    on_circle = np.exp(1j * np.radians(data))
    mean = np.mean(on_circle)
    mean_angle = np.angle(mean, deg=True)
    around_zero = np.angle(on_circle / mean, deg=True)
    min, max = np.percentile(around_zero, percentiles[0], method="lower"), np.percentile(around_zero, percentiles[1], method="higher")
    return float(angular_modulo(min + mean_angle)), float(angular_modulo(max + mean_angle))

def angular_modular_kde(data:Union[pl.Series, np.ndarray], bw_factor = 1) -> scipy.stats.gaussian_kde:
    """
    Computes KDE with 180 degrees repeated padding on both sides - range is -360...360, so the range -180...180 should have somewhat reliable likelihoods
    """
    if isinstance(data, pl.Series):
        data = data.drop_nulls().to_numpy()

    data = np.sort(data)
    data_left = data[data < 0] + 360
    data_right = data[data > 0] - 360
    return scipy.stats.gaussian_kde(np.concatenate([ data_left, data, data_right ]))

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
    if not (m := re.match("^hb_(\\d+)_", col)):
        return None
    ix = int(m.group(1))
    if len(hbonds) <= ix: return None

    # cWB has two identical bonds, because the nitrogen has two hydrogens...
    if hbonds.count(hbonds[ix]) > 1:
        disambig = " " + chr(ord('α') + hbonds[:ix].count(hbonds[ix]))
    else:
        disambig = ""

    if col.endswith("_donor_angle"):
        if len(hbonds) <= ix: return None
        return format_angle_label(hbonds[ix][:3], swap=swap)+disambig
    elif col.endswith("_acceptor_angle"):
        if len(hbonds) <= ix: return None
        return format_angle_label(hbonds[ix][1:], swap=swap)+disambig
    else:
        return format_length_label(hbonds[ix][1:3], swap=swap)+disambig

def is_swapped(pair_type: PairType):
    symtype = pair_type.type.lower() in ["cww", "tww", "chh", "thh", "css", "tss"]
    return pair_type.bases_str == "C-G" and symtype or pair_type.type.lower() not in pair_defs.pair_families
def format_pair_type(pair_type: PairType, is_dna = False, is_rna=False):
    pair_kind, pair_bases = pair_type.to_tuple()
    if is_dna == True:
        pair_bases = pair_bases.replace("U", "T")
    elif is_rna == True:
        pair_bases = pair_bases.replace("T", "U")
    if is_swapped(pair_type):
        assert not is_swapped(pair_type.swap())
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

def draw_pair_img_highligh(ax, img, highlight: Optional[pl.DataFrame]):
    if img is None:
        return
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
        # fig.suptitle(title + f" ({len(df)}, class {determine_bp_class(df, pair_type)})")
        fig.suptitle(title + f" ({len(df)} observations)")

    plot_offset = 0 if images is None else 1
    for i, h in enumerate(hs):
        make_histogram_group(dataframes, [ p[1][i+plot_offset] for p in pages ], [h.title] * len(dataframes), pair_type, h)
    if images is not None:
        for p, img, highlight in zip(pages, images, highlights or itertools.repeat(None)):
            ax = p[1][0]
            draw_pair_img_highligh(ax, img, highlight)
    if highlights is not None:
        for p, highlight in zip(pages, highlights):
            if highlight is None or len(highlight) == 0:
                continue

            for ax, h in zip(p[1][plot_offset:], hs):
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

def make_pairplot_page(df: pl.DataFrame, outdir: str, pair_type: PairType, variables: list[pl.Expr], title_suffix = ""):
    title = f"{format_pair_type(pair_type)} Pairplot{title_suffix}"
    grid: sns.PairGrid = sns.pairplot(df.select(*variables).to_pandas(),
                                        kind="scatter" if len(df) < 400 else "kde",
                                        #  kind="hist"
                                        )
    grid.figure.suptitle(title)

    yield save(grid.figure, title, outdir)

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
            pl.col("^hb_\\d+_donor_angle$").degrees(),
            pl.col("^hb_\\d+_acceptor_angle$").degrees(),
        )
    if "C1_C1_euler_phi" in df.columns:
        df = df.with_columns(
            C1_C1_euler_phicospsi=
                (pl.col("C1_C1_euler_phi") - pl.col("C1_C1_euler_psi") * pl.col("C1_C1_euler_theta").radians().cos() + 180 ) % 360 - 180,
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

def get_kde_mode(kde: scipy.stats.gaussian_kde, data: np.ndarray):
    data = np.sort(data)
    range = np.min(data), np.max(data)
    coarse_space = np.linspace(range[0], range[1], num=200)
    coarse_argmax = np.argmax(kde.logpdf(coarse_space))
    coarse_argmax = min(len(coarse_space)-2, max(1, coarse_argmax))

    fine_range = coarse_space[coarse_argmax-1], coarse_space[coarse_argmax+1]
    fine_data = data[(data <= fine_range[1]) & (data >= fine_range[0])]
    if len(fine_data) == 0:
        return coarse_space[coarse_argmax]
    data_argmax = np.argmax(kde.logpdf(fine_data))
    return float(fine_data[data_argmax])

def calculate_stats(df: pl.DataFrame, pair_type, skip_kde: bool):
    if len(df) == 0:
        raise ValueError("No data")
    def make_kde(data, cname):
        if is_angular_modular(cname):
            kde = angular_modular_kde(sample_for_kde(data))
        else:
            kde = scipy.stats.gaussian_kde(sample_for_kde(data))
        kde.set_bandwidth(kde.scotts_factor() * 1.5)
        return kde
    def calc_mean(data, cname):
        if is_angular_modular(cname):
            return angular_modular_mean(sample_for_kde(data))
        else:
            return float(np.mean(data))
    def calc_std(data, cname, mean=None):
        if is_angular_modular(cname):
            return angular_modular_std(sample_for_kde(data), mean=mean)
        else:
            return float(np.std(data))

    def calc_diff(data, cname):
        if is_angular_modular(cname):
            return angular_modulo(data)
        else:
            return data

    # columns = df.select(pl.col("^hb_\\d+_(length|donor_angle|acceptor_angle)$")).columns
    # columns = df.select(pl.col("^hb_\\d+_length$"), pl.col("^C1_C1_(yaw|pitch|roll)(1|2)$")).columns
    columns = df.select(pl.col("^hb_\\d+_(length|donor_angle|acceptor_angle)$"), pl.col("^C1_C1_(yaw|pitch|roll)(1|2)$")).columns
    # print(f"{columns=}")
    cdata = [ df[c].drop_nulls().to_numpy() for c in columns ]
    medians = [ float(np.median(c)) if len(c) > 0 and not is_angular_modular(cname) else None for c, cname in zip(cdata, columns) ]
    means = [ calc_mean(c, cname) if len(c) > 0 else None for c, cname in zip(cdata, columns) ]
    stds = [ calc_std(c, cname) if len(c) > 1 else None for c, cname in zip(cdata, columns) ]

    if skip_kde:
        kdes = [ None ] * len(columns)
    else:
        kdes = [ make_kde(c, cname) if len(c) > 5 else None for c, cname in zip(cdata, columns) ]

    kde_modes = [ None if kde is None else float(get_kde_mode(kde, c)) for kde, c in zip(kdes, cdata) ]
    kde_mode_stds = [ None if kde_mode is None else calc_std(c, cname, mean=kde_mode) for kde_mode, c, cname in zip(kde_modes, cdata, columns) ]
    def calc_datapoint_mode_deviations(df):
        ms = [ median if kde_mode is None else kde_mode for kde_mode, median in zip(kde_modes, medians) ]
        ss = [ std if kde_std is None else kde_std for kde_std, std in zip(kde_mode_stds, stds) ]
        return np.sum([
            (calc_diff(df[c] - kde_mode, c).abs() / abs(kde_mode_std)).fill_null(10).to_numpy()
            for kde_mode, kde_mode_std, c in zip(ms, ss, columns)
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
        "nicest_bp": list(next(df[nicest_basepair, ["pdbid", "model", "chain1", "res1", "nr1", "ins1", "alt1", "chain2", "res2", "nr2", "ins2", "alt2"]].iter_rows())),
        "nicest_bp_index": nicest_basepair,
        "nicest_bp_indices": nicest_basepairs,
        # **{
        #     f"hb_{i}_label": get_label(f"hb_{i}_length", pair_type) for i in range(len(pair_defs.get_hbonds(pair_type)))
        # },
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
    p = subprocess.run(command, capture_output=True, timeout=30)
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

def reexport_df(df: pl.DataFrame, columns, drop: list[str] = []):
    df = df.with_columns(
        is_some_quality.alias("jirka_approves"),
        *[
            pl.Series(columns[col](df)).alias(col)
            for col in columns
        ]
    )
    df = df.drop([col for col in df.columns if re.match(r"[DR]NA-(0-1[.]8|1[.]8-3[.]5)(-r\d+)?", col)])
    df = df.drop([col for col in df.columns if col.startswith('^') in drop])
    df = df.drop(["label"])
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
            print(f"{file}: {len(df)} rows, types: {dict(sorted([ (str(pair_defs.PairType.from_tuple(k)), len(gdf)) for k, gdf in groups ], key=lambda x: x[1], reverse=True))}")
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
    pt_family_dict = { pt: ix + 1 for ix, pt in enumerate(pair_defs.pair_families) }
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
        df2.write_excel(workbook, worksheet="All", dtype_formats={ pl.Float64: "0.00", pl.Float32: "0.00" }, hidden_columns=["hb_ix", "Pair bases"])
        resolutions = df2["Resolution cutoff"].unique().to_list()
        for r in resolutions:
            df2.filter(pl.col("Resolution cutoff") == r)\
                .write_excel(workbook, worksheet=f"{r}", dtype_formats={ pl.Float64: "0.00", pl.Float32: "0.00" }, hidden_columns=["hb_ix", "Pair bases", "Resolution cutoff"])


def round_boundary(column: str, down: bool, value: Optional[float]):
    if value is None:
        return None
    value = value if down else -value
    if column.endswith("_length") or column.startswith("rmsd_") or column in ["coplanarity_shift1", "coplanarity_shift2"]:
        value = math.floor(value * 10) / 10
    else:
        value = math.floor(value / 5) * 5
    return value if down else -value

def calculate_boundaries(df: pl.DataFrame, pair_type: PairType):
    hbonds = pair_defs.get_hbonds(pair_type)
    hb_is_o2prime = [ "O2'" in hb for hb in hbonds ]
    blacklisted_columns = [
        *[ f"hb_{i}_OOPA1" for i in range(len(hbonds)) if hb_is_o2prime[i] ],
        *[ f"hb_{i}_OOPA2" for i in range(len(hbonds)) if hb_is_o2prime[i] ],
    ]
    boundary_columns = {
        "yaw1": "C1_C1_yaw1",
        "pitch1": "C1_C1_pitch1",
        "roll1": "C1_C1_roll1",
        "yaw2": "C1_C1_yaw2",
        "pitch2": "C1_C1_pitch2",
        "roll2": "C1_C1_roll2",
        "coplanarity_angle": "coplanarity_angle",
        "coplanarity_edge_angle1": "coplanarity_edge_angle1",
        "coplanarity_edge_angle2": "coplanarity_edge_angle2",
        "coplanarity_shift1": "coplanarity_shift1",
        "coplanarity_shift2": "coplanarity_shift2",
        "hb_0_acceptor_angle": "hb_0_acceptor_angle",
        "hb_1_acceptor_angle": "hb_1_acceptor_angle",
        "hb_2_acceptor_angle": "hb_2_acceptor_angle",
        "hb_0_donor_angle": "hb_0_donor_angle",
        "hb_1_donor_angle": "hb_1_donor_angle",
        "hb_2_donor_angle": "hb_2_donor_angle",
        "hb_0_OOPA1": "hb_0_OOPA1",
        "hb_0_OOPA2": "hb_0_OOPA2",
        "hb_1_OOPA1": "hb_1_OOPA1",
        "hb_1_OOPA2": "hb_1_OOPA2",
        "hb_2_OOPA1": "hb_2_OOPA1",
        "hb_2_OOPA2": "hb_2_OOPA2",
    }
    pt_family_dict = { pt: ix + 1 for ix, pt in enumerate(pair_defs.pair_families) }
    def calc_boundary(col, df: pl.DataFrame=df):
        # return df[col].min(), df[col].max()
        if col in blacklisted_columns:
            min, max = None, None
        elif len(df) <= 1 or df[col].is_null().mean() > 0.1:
            min, max = None, None
        
        elif is_angular_modular(col):
            min, max = angular_modular_minmax(df[col].drop_nulls()) or [None, None]
            # min, max = angular_modular_minmax(df[col].drop_nulls(), percentiles=(1, 99)) or [None, None]
        else:
            min, max = df[col].min(), df[col].max()
            # min, max = df[col].quantile(0.001, "lower"), df[col].quantile(0.999, "higher")
        
        min = round_boundary(col, True, min)
        max = round_boundary(col, False, max)

        # return pl.Series([], dtype=pl.Float32)
        return pl.Series([min, max], dtype=pl.Float32)
    
    # hb_dict = dict(enumerate(hbonds))

    max_lengths = [
        4.2 if pair_defs.is_ch_bond(pair_type, hb) else 4.0
        for hb in hbonds
    ]
    max_lengths = dict(enumerate(max_lengths))

    boundaries = pl.DataFrame({
        "family_id": [ pt_family_dict.get(pair_type.type.lower(), 99) ] * 2,
        "family": [ pair_type.full_family ] * 2,
        "bases": [ pair_type.bases_str ] * 2,
        "count": [ len(df) ] * 2,
        "boundary": [ "min", "max" ],
        **{
            key: calc_boundary(col)
            for key, col in boundary_columns.items()
        },
        "hb_0_length": [ None, max_lengths.get(0) ],
        "hb_1_length": [ None, max_lengths.get(1) ],
        "hb_2_length": [ None, max_lengths.get(2) ],
        "min_bond_length": [ None, 3.8 ],
    })
    return boundaries

def main(argv):
    import argparse
    parser = argparse.ArgumentParser(description="Generate histogram plots")
    parser.add_argument("input_file", help="Input file", nargs="+")
    parser.add_argument("--residue-directory", help="Directory with residue lists, used to select representative set residues. Currently, lists RNA-1.8-3.5, RNA-0-1.8, DNA-1.8-3.5, DNA-0-1.8 are expected.")
    parser.add_argument("--reexport", default='none', choices=['none', 'partitioned'], help="Write out parquet files with calculated statistics columns (log likelihood, mode deviations)")
    parser.add_argument("--include-nears", default=False, action="store_true", help="If FR3D is run in basepair_detailed mode, it reports near basepairs (denoted as ncWW). By default, we ignore them, but this option includes them in the output.")
    parser.add_argument("--filter-pair-type", default=None, help="Comma separated list of pair types to include in the result (formatted as cWW-AC). By default all are included.")
    parser.add_argument("--skip-kde", default=False, action="store_true", help="Skip generating kernel density estimates for histograms, image selection and 'niceness' score calculation.")
    parser.add_argument("--skip-image", default=False, action="store_true", help="Skip generating images for the nicest basepairs (use if the gen_contact_images.py script is broken)")
    parser.add_argument("--drop-columns", default=[], nargs="*", help="remove the specified columns from the reexport output (regex supported when in ^...$)")
    parser.add_argument("--output-dir", "-o", required=True, help="Output directory")
    args = parser.parse_args(argv)

    only_pairtypes = None if args.filter_pair_type is None else set(pair_defs.PairType.parse(pt) for pt in args.filter_pair_type.split(","))

    if args.residue_directory:
        residue_lists = residue_filter.read_id_lists(args.residue_directory)
    else:
        residue_lists = None
        global is_some_quality # TODO: configurable resolution cutoffs
        is_some_quality = pl.lit(True)
    
    if args.skip_kde:
        global hist_kde
        hist_kde = False

    results = []
    boundaries = []
    all_statistics = []

    os.makedirs(args.output_dir, exist_ok=True)
    for pair_type, df in enumerate_pair_types(args.input_file, args.include_nears):
        if only_pairtypes is not None and pair_type.without_n() not in only_pairtypes:
            print(f"Skipping {pair_type} because it is not in {only_pairtypes}")
            continue
        
        if residue_lists:
            df = residue_filter.add_res_filter_columns(df, residue_lists)
        else:
            df = df.with_columns(**{
                "RNA-0-1.8": is_rna & (pl.col("resolution") <= 1.8),
                "RNA-1.8-3.5": is_rna & (pl.col("resolution") <= 3.5) & (pl.col("resolution") > 1.8),
                "DNA-0-1.8": is_dna & (pl.col("resolution") <= 1.8),
                "DNA-1.8-3.5": is_dna & (pl.col("resolution") <= 3.5) & (pl.col("resolution") > 1.8),
            })
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
            stat_columns, stats = calculate_stats(dff, pair_type, args.skip_kde)
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
            if args.skip_image:
                basepair_images = None
            else:
                try:
                    basepair_images = [ create_pair_image(dff[bp], args.output_dir, pair_type) if bp >= 0 else None for bp in nicest_bps ] * len(resolutions) if nicest_bps is not None else []
                except Exception as e:
                    print(f"Error generating images for {pair_type}: {e}")
                    print(f"If the image generation doesn't work (e.g. PyMOL is not installed), you can skip it with --skip-image flag")
                    exit(1)
            dna_rna_highlights = [ dff[bp] if bp >= 0 else None for bp in nicest_bps ] if nicest_bps is not None else []
            output_files = []
            output_files.extend(
                f for f in make_bond_pages(df, args.output_dir, pair_type, hbond_histogram_defs, images=basepair_images, highlights=dna_rna_highlights, title_suffix=" - H-bonds"
                )
            )
            output_files.extend(
                make_bond_pages(dff, args.output_dir, pair_type, coplanarity_histogram_defs, highlights=dna_rna_highlights, title_suffix= " - Coplanarity")
            )
            output_files.extend(
                make_bond_pages(dff, args.output_dir, pair_type, coplanarity_histogram_defs2, highlights=dna_rna_highlights, title_suffix= " - coplanarity2")
            )
            output_files.extend(
                make_bond_pages(dff, args.output_dir, pair_type, rmsd_histogram_defs, highlights=dna_rna_highlights, title_suffix= " - RMSD to nicest BP")
            )
            # output_files.extend(
            #     make_pairplot_page(dff, args.output_dir, pair_type, variables=[
            #         pl.col(f"C1_C1_euler_phi").alias("Euler Φ"),
            #         pl.col(f"C1_C1_euler_theta").alias("Euler Θ"),
            #         pl.col(f"C1_C1_euler_psi").alias("Euler Ψ"),
            #         pl.col(f"C1_C1_euler_phicospsi").alias("Euler Φ-cos(Θ)Ψ"),
            #         pl.col(f"C1_C1_yaw1").alias("Yaw 1"),
            #         pl.col(f"C1_C1_pitch1").alias("Pitch 1"),
            #         pl.col(f"C1_C1_roll1").alias("Roll 1"),
            #         pl.col(f"C1_C1_yaw2").alias("Yaw 2"),
            #         pl.col(f"C1_C1_pitch2").alias("Pitch 2"),
            #         pl.col(f"C1_C1_roll2").alias("Roll 2"),
            #     ], title_suffix=" - Various N1-C1' reference frame angles")
            # )
            # output_files.extend(
            #     make_pairplot_page(dff, args.output_dir, pair_type, variables=[
            #         pl.col(f"coplanarity_angle").alias("Plane normal angle"),
            #         pl.col(f"C1_C1_yaw1").alias("Yaw 1"),
            #         pl.col(f"C1_C1_pitch1").alias("Pitch 1"),
            #         pl.col(f"C1_C1_roll1").alias("Roll 1"),
            #         pl.col(f"coplanarity_edge_angle1").alias("Edge 1 / plane 2"),
            #         pl.col(f"coplanarity_edge_angle2").alias("Edge 2 / plane 1"),
            #         pl.col("hb_1_OOPA1").alias("H-bond 2 / plane 1"),
            #         pl.col("hb_1_OOPA2").alias("H-bond 2 / plane 2"),
            #         pl.col("coplanarity_shift1").alias("Edge 1/plane 2 shift"),
            #         pl.col("coplanarity_shift2").alias("Edge 2/plane 1 shift"),

            #     ], title_suffix=" - Other coplanarity metrics")
            # )
            # output_files.extend(
            #     make_pairplot_page(dff, args.output_dir, pair_type, variables=[
            #         pl.col(f"C1_C1_yaw1").alias("Yaw 1"),
            #         pl.col(f"C1_C1_pitch1").alias("Pitch 1"),
            #         pl.col(f"C1_C1_roll1").alias("Roll 1"),
            #         pl.col(f"rmsd_C1N_frames"),
            #         pl.col(f"rmsd_edge1"),
            #         pl.col(f"rmsd_edge2"),
            #         pl.col(f"rmsd_edge_C1N_frame"),
            #         pl.col(f"rmsd_all_base"),
            #     ], title_suffix=" - Various N1-C1' reference frame angles")
            # )
            # output_files.extend(
            #     f
            #     for column in [0, 1, 2]
            #     for f in make_bond_pages(df, args.output_dir, pair_type, [ h.select_columns(column) for h in hbond_histogram_defs], images=dna_rna_images, highlights=dna_rna_highlights, title_suffix=f" #{column}")
            # )
            all_statistics.extend(statistics)
            hb_filters = [
                # pl.col(f"hb_{i}_length") <= 4.0 if "C" in hb[1] or "C" in hb[2] else pl.col(f"hb_{i}_length") <= 3.8
                # for i, hb in enumerate(pair_defs.get_hbonds(pair_type, throw=False))
                # if f"hb_{i}_length" in df.columns
            ]
            boundaries.append(calculate_boundaries(
                df.filter(
                    pl.all_horizontal(is_some_quality, *hb_filters)
                ),
                pair_type
            ))
        if args.reexport == "partitioned":
            reexport_df(df, stat_columns or [], drop=args.drop_columns).write_parquet(os.path.join(args.output_dir, f"{pair_type.normalize_capitalization()}.parquet"))
            reexport_df(df.filter(is_some_quality), stat_columns or [], drop=args.drop_columns).write_parquet(os.path.join(args.output_dir, f"{pair_type.normalize_capitalization()}-filtered.parquet"))
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
    # results.sort(key=lambda r: (r["bp_class"] or 5, pair_defs.PairType.from_tuple(r["pair_type"])))
    results.sort(key=lambda r: (1, pair_defs.PairType.from_tuple(r["pair_type"])))
    output_files = [ f for r in results for f in r["files"] ]

    subprocess.run(["gs", "-dBATCH", "-dNOPAUSE", "-q", "-sDEVICE=pdfwrite", "-dPDFSETTINGS=/prepress", f"-sOutputFile={os.path.join(args.output_dir, 'plots-merged.pdf')}", *output_files])
    print("Wrote", os.path.join(args.output_dir, 'plots-merged.pdf'))
    boundaries_df: pl.DataFrame = pl.concat(boundaries)
    boundaries_df = boundaries_df.sort("family_id", "bases", "family", "boundary", descending=[False, False, False, True])
    boundaries_df.write_csv(os.path.join(args.output_dir, "boundaries.csv"))

    boundaries_reformat = boundaries_df.group_by("family_id", "family", "bases").agg(*itertools.chain(*[
        [pl.col(col).filter(pl.col("boundary") == "min").first().alias(col + "_min"), pl.col(col).filter(pl.col("boundary") == "max").first().alias(col + "_max")]
        for col in boundaries_df.columns
        if col not in ["family", "family_id", "bases", "count", "boundary"]
    ])).sort("family_id", "bases", "family")
    boundaries_reformat.write_csv(os.path.join(args.output_dir, "boundaries2.csv"))
    import xlsxwriter
    with xlsxwriter.Workbook(os.path.join(args.output_dir, "boundaries.xlsx")) as workbook:
        boundaries_df.write_excel(workbook, worksheet="Narrow long table", dtype_formats={ pl.Float64: "0.00", pl.Float32: "0.00" })
        boundaries_reformat.write_excel(workbook, worksheet="Wider shorter table", dtype_formats={ pl.Float64: "0.00", pl.Float32: "0.00" })
    # save_statistics(all_statistics, args.output_dir)

    with open(os.path.join(args.output_dir, "output.json"), "w") as f:
        import json
        json.dump(results, f, indent=4)

if __name__ == "__main__":
    main(sys.argv[1:])
