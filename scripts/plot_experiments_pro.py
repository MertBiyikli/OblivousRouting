#!/usr/bin/env python3
from __future__ import annotations

import re
from collections.abc import Callable
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker


# ======================
# CONFIG
# ======================
RESULTS_CSV = Path("/Users/halilibrahim/Desktop/Thesis/ObliviousRouting/results/combined.csv")
OUT_DIR = Path("/Users/halilibrahim/Desktop/Thesis/ObliviousRouting/results/2dgrids/")
OUT_DIR.mkdir(parents=True, exist_ok=True)

ORACLE_TIME_COL = "avg_oracle_time_micro_seconds"

TIME_COLUMN_ALIASES: dict[str, tuple[str, float]] = {
    "total_time_micro_seconds": ("total_time_ms", 1_000.0),
    "solve_time_micro_seconds": ("solve_time_ms", 1_000.0),
    "transformation_time_micro_seconds": ("transformation_time_ms", 1_000.0),
    "avg_oracle_time_micro_seconds": ("avg_oracle_time_ms", 1_000.0),
    "mendel_total_micro_seconds": ("mendel_total_ms", 1_000.0),
    "mendel_avg_micro_seconds": ("mendel_avg_ms", 1_000.0),
}

# Typical LaTeX paper sizes
FIGSIZE_SINGLE = (2.87, 2.17)   # thesis single-column-ish
FIGSIZE_DOUBLE = (5.91, 2.17)     # thesis double-column-ish

EXPORT_PNG = True
DPI_RASTER = 450

SHOW_ERROR_BARS = False
LOG_EPS = 1e-3  # microseconds, safe lower bound for log-plots

# Legend control
LEGEND_NCOL = 2
LEGEND_OUTSIDE = True  # put legend above plot to avoid occluding data
# Set False to strip legends from every individual plot and instead emit a
# single standalone legend_plate.pdf/png that can be included once in a paper.
EMBED_LEGEND = False
# Set False to strip all axis labels, titles, and named tick labels from every
# plot.  Axis/scale tick numbers are kept; all human-readable text is removed
# so plots can be annotated via LaTeX captions / pgfplots overlay instead.
EMBED_AXIS_LABELS = False


# ======================
# STYLE (publication-grade)
# ======================

def set_paper_style():
    mpl.rcParams.update({
        # ── Typography ──────────────────────────────────────────────────────────
        "font.family": "serif",
        "font.serif": ["Times New Roman", "Times", "STIXGeneral", "DejaVu Serif"],
        "mathtext.fontset": "stix",
        "font.size": 9,
        "axes.labelsize": 8,
        "axes.titlesize": 8,
        "axes.titlepad": 4,
        "legend.title_fontsize": 7,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,

        # ── Lines / markers ──────────────────────────────────────────────────
        "lines.linewidth": 0.5,
        "lines.markersize": 2.5,
        "lines.markeredgewidth": 0.4,

        # ── Axes & spines ────────────────────────────────────────────────────
        "axes.linewidth": 0.8,
        "axes.labelpad": 3.0,
        "axes.spines.top": False,
        "axes.spines.right": False,

        # ── Ticks ────────────────────────────────────────────────────────────
        "xtick.major.width": 0.7,
        "ytick.major.width": 0.7,
        "xtick.minor.width": 0.5,
        "ytick.minor.width": 0.5,
        "xtick.major.size": 3.0,
        "ytick.major.size": 3.0,
        "xtick.minor.size": 1.8,
        "ytick.minor.size": 1.8,
        "xtick.major.pad": 2.5,
        "ytick.major.pad": 2.5,
        "xtick.direction": "in",
        "ytick.direction": "in",

        # ── Grid ─────────────────────────────────────────────────────────────
        "axes.grid": True,
        "axes.axisbelow": True,        # grid behind data
        "grid.color": "#CCCCCC",
        "grid.alpha": 0.6,
        "grid.linewidth": 0.4,
        "grid.linestyle": "--",

        # ── Legend ───────────────────────────────────────────────────────────
        "legend.frameon": False,
        "legend.borderpad": 0.3,
        "legend.labelspacing": 0.3,
        "legend.handlelength": 1.6,
        "legend.handleheight": 0.7,
        "legend.handletextpad": 0.4,
        "legend.columnspacing": 0.8,
        "legend.fontsize": 7,

        # ── Figure & output ──────────────────────────────────────────────────
        "figure.dpi": 150,
        "pdf.fonttype": 42,            # embed TrueType (required by IEEE/ACM)
        "ps.fonttype": 42,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.03,
        "savefig.dpi": 450,
    })


def savefig_all(fig: plt.Figure, outbase: Path):
    fig.savefig(outbase.with_suffix(".pdf"))
    if EXPORT_PNG:
        fig.savefig(outbase.with_suffix(".png"), dpi=DPI_RASTER)


def normalize_time_columns_to_microseconds(df: pd.DataFrame) -> pd.DataFrame:
    """
    Accept both legacy millisecond CSVs and the newer microsecond CSVs.

    For each canonical ``*_micro_seconds`` column:
      * keep the microsecond values when already present,
      * otherwise backfill from the corresponding ``*_ms`` column × 1000,
      * create the canonical column when only the legacy millisecond column exists.
    """
    df = df.copy()

    for canonical_col, (legacy_col, scale) in TIME_COLUMN_ALIASES.items():
        if canonical_col in df.columns:
            df[canonical_col] = pd.to_numeric(df[canonical_col], errors="coerce")
            if legacy_col in df.columns:
                legacy_vals = pd.to_numeric(df[legacy_col], errors="coerce") * scale
                df[canonical_col] = df[canonical_col].fillna(legacy_vals)
        elif legacy_col in df.columns:
            df[canonical_col] = pd.to_numeric(df[legacy_col], errors="coerce") * scale

    return df


# ======================
# DATA AGGREGATION
# ======================

def aggregate_mean_std(df: pd.DataFrame, ycol: str) -> pd.DataFrame:
    g = (
        df.groupby(["solver", "num_edges"], as_index=False)[ycol]
        .agg(["mean", "std", "count"])
        .reset_index()
    )
    g.rename(columns={"count": "n"}, inplace=True)
    g["std"] = g["std"].fillna(0.0)
    return g


# ======================
# SOLVER LABELS + ORDER
# ======================

def pretty_solver_name(s: str) -> str:
    # If it already has a nice format with parentheses, return as-is
    if "(" in s and ")" in s:
        return s

    # Otherwise, apply the old mapping for legacy format
    mapping = {
        "electrical": "Electrical Flow",
        "electrical_parallel": "Electrical Flow (parallel)",
        "raecke_frt": "Räcke-FRT",
        "raecke_ckr": "Räcke–Fast-CKR",
        "raecke_mst": "Räcke–MST",
        "cohen": "LP (Applegate–Cohen)",
        "lp": "LP (Applegate–Cohen)",
        "random_mst": "Räcke–MST",
        "mst": "Räcke–MST",
        "frt": "Räcke-FRT",
        "ckr": "Räcke-Fast-CKR",
        "frt_mendel": "Räcke-FRT (Mendel)",
        "ckr_mendel": "Räcke-Fast-CKR (Mendel)",
        "raecke_frt_mendel": "Räcke-FRT (Mendel)",
        "raecke_ckr_mendel": "Räcke–Fast-CKR (Mendel)",
        "raecke_random_mst": "Räcke–MST",
    }
    return mapping.get(s, s)


def solver_sort_key(s: str) -> tuple:
    """
    Global solver order used across *all* plots for consistency.
    electrical → CKR → CKR (Mendel) → FRT → FRT (Mendel) → MST → ... → LP (always last)
    Handles both old format (electrical) and new format with variants (Electrical Flow (naive))
    """
    # Normalize solver name - extract base name if it has variants in parentheses
    base_s = s.split("(")[0].strip().lower().replace(" ", "_")

    # Handle new format names
    if "electrical" in base_s:
        if "parallel" in s.lower():
            order_val = 1
        else:
            order_val = 0  # Both naive and sketching go first
    elif "ckr" in base_s:
        if "mendel" in s.lower():
            order_val = 3
        else:
            order_val = 2
    elif "frt" in base_s:
        if "mendel" in s.lower():
            order_val = 5
        else:
            order_val = 4
    elif any(x in base_s for x in ["mst", "random"]):
        order_val = 6
    elif any(x in base_s for x in ["cohen", "lp", "applegate"]):
        order_val = 99
    else:
        order_val = 90

    return (order_val, s)


# ======================
# COLORS (unique per solver, consistent across plots)
# ======================

# Curated colorblind-safe palette (Wong 2011 + Tol extensions).
# Stays readable in grayscale and on projectors.
_PALETTE = [
    "#0072B2",  # blue
    "#D55E00",  # vermillion
    "#009E73",  # green
    "#CC79A7",  # pink
    "#E69F00",  # orange
    "#56B4E9",  # sky-blue
    "#F0E442",  # yellow
    "#000000",  # black
    "#999999",  # grey
    "#882255",  # wine
    "#117733",  # dark green
    "#44AA99",  # cyan
    "#DDCC77",  # khaki
    "#AA4499",  # purple
]

# Fixed color mapping for all possible solvers to ensure consistency across all runs
_SOLVER_COLOR_MAP = {
    # Electrical Flow variants
    "Electrical Flow (naive)": "#0072B2",           # blue
    "Electrical Flow (sketching)": "#D55E00",       # vermillion
    "Electrical Flow (parallel)": "#E69F00",        # orange
    
    # Raecke CKR - Flat HST
    "Raecke CKR (Flat HST)": "#009E73",             # green
    "Raecke CKR + MendelScaling (Flat HST)": "#56B4E9",  # sky-blue
    
    # Raecke FRT - Flat HST
    "Raecke FRT (Flat HST)": "#CC79A7",             # pink
    "Raecke FRT + MendelScaling (Flat HST)": "#F0E442",  # yellow
    
    # Random MST - Flat HST
    "Random MST (Flat HST)": "#000000",             # black
    
    # Raecke CKR - Pointer HST
    "Raecke CKR (Pointer HST)": "#999999",          # grey
    "Raecke CKR + MendelScaling (Pointer HST)": "#882255",  # wine
    
    # Raecke FRT - Pointer HST
    "Raecke FRT (Pointer HST)": "#117733",          # dark green
    "Raecke FRT + MendelScaling (Pointer HST)": "#44AA99",  # cyan
    
    # Random MST - Pointer HST
    "Random MST (Pointer HST)": "#DDCC77",          # khaki
    
    # LP
    "LP Applegate-Cohen": "#AA4499",                # purple
}

_LINESTYLES = ["-", "--", "-.", ":", (0, (3, 1, 1, 1)), (0, (5, 1))]


def solver_palette_unique(solvers: list[str]) -> dict[str, str]:
    """
    Unique color per solver using predefined mapping for consistency.
    Falls back to dynamic palette if solver not in mapping.
    """
    colors = {}
    for i, s in enumerate(solvers):
        if s in _SOLVER_COLOR_MAP:
            colors[s] = _SOLVER_COLOR_MAP[s]
        else:
            # Fallback: use dynamic palette for unknown solvers
            colors[s] = _PALETTE[i % len(_PALETTE)]
    return colors


def solver_linestyles(solvers: list[str]) -> dict[str, str]:
    """Unique linestyle per solver (helps grayscale / print)."""
    return {s: _LINESTYLES[i % len(_LINESTYLES)] for i, s in enumerate(solvers)}


def solver_markers(solvers: list[str]) -> dict[str, str]:
    """
    Optional: different markers help when printed in grayscale.
    Still keeps color unique as requested.
    """
    pool = ["o", "s", "^", "D", "v", "P", "X", ">", "<", "h", "*"]
    return {s: pool[i % len(pool)] for i, s in enumerate(solvers)}


# ======================
# PLOTTING
# ======================

def add_legend(ax: plt.Axes, handles, labels):
    if not labels or not EMBED_LEGEND:
        return

    if LEGEND_OUTSIDE:
        ax.legend(
            handles, labels,
            loc="lower left",
            bbox_to_anchor=(0.0, 1.01),
            ncol=min(LEGEND_NCOL, max(1, len(labels))),
            frameon=False,
            borderaxespad=0.0,
            columnspacing=0.8,
            handletextpad=0.5,
        )
    else:
        ax.legend(frameon=False, ncol=min(LEGEND_NCOL, max(1, len(labels))))


def _apply_labels(ax: plt.Axes,
                  xlabel: str = "",
                  ylabel: str = "",
                  title:  str | None = None):
    """Set axis / title labels only when EMBED_AXIS_LABELS is True."""
    if not EMBED_AXIS_LABELS:
        return
    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title)


def plot_lines(
        df_agg: pd.DataFrame,
        solvers: list[str],
        colors: dict[str, str],
        markers: dict[str, str],
        xcol: str,
        y_mean_col: str,
        y_std_col: str,
        xlabel: str,
        ylabel: str,
        xlog: bool,
        ylog: bool,
        figsize,
        outpath: Path,
        linestyles: dict[str, str] | None = None,
):
    fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)

    handles = []
    labels = []

    for s in solvers:
        sub = df_agg[df_agg["solver"] == s].sort_values(xcol)
        if sub.empty:
            continue

        x = sub[xcol].to_numpy()
        y = sub[y_mean_col].to_numpy()
        e = sub[y_std_col].to_numpy()

        if ylog:
            y = np.maximum(y, LOG_EPS)
            e = np.minimum(e, y - LOG_EPS)

        ls = linestyles[s] if linestyles else "-"
        h, = ax.plot(
            x, y,
            linestyle=ls,
            marker=markers[s],
            markersize=2.5,
            markeredgewidth=0.4,
            markeredgecolor="white",
            label=pretty_solver_name(s),
            color=colors[s],
        )
        handles.append(h)
        labels.append(pretty_solver_name(s))

        if SHOW_ERROR_BARS and np.any(e > 0):
            ax.fill_between(x, y - e, y + e, alpha=0.10, color=colors[s], linewidth=0)

    if xlog:
        ax.set_xscale("log")
    if ylog:
        ax.set_yscale("log")

    _apply_labels(ax, xlabel, ylabel)
    ax.minorticks_on()

    add_legend(ax, handles, labels)

    savefig_all(fig, outpath)
    plt.close(fig)

def aggregate_mean_std_by_solver(df: pd.DataFrame, ycol: str) -> pd.DataFrame:
    g = (
        df.groupby(["solver"], as_index=False)[ycol]
        .agg(["mean", "std", "count"])
        .reset_index()
    )
    g.rename(columns={"count": "n"}, inplace=True)
    g["std"] = g["std"].fillna(0.0)
    return g


def plot_solver_average_bars(
        df_solver_agg: pd.DataFrame,
        solvers: list[str],
        colors: dict[str, str],
        y_mean_col: str,
        y_std_col: str,
        ylabel: str,
        ylog: bool,
        figsize,
        outpath: Path,
):
    fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)

    sub = df_solver_agg.set_index("solver").reindex(solvers).reset_index()

    y = sub[y_mean_col].to_numpy(dtype=float)
    e = sub[y_std_col].to_numpy(dtype=float)

    if ylog:
        y = np.maximum(y, LOG_EPS)
        e = np.minimum(e, y - LOG_EPS)

    x = np.arange(len(solvers))
    ax.bar(
        x, y,
        color=[colors[s] for s in solvers],
        edgecolor="white",
        linewidth=0.5,
        alpha=0.88,
        zorder=3,
    )

    # Always show error bars for bars (they carry statistical meaning)
    ax.errorbar(
        x, y, yerr=e,
        fmt="none", ecolor="#333333",
        capsize=2.5, capthick=0.7,
        elinewidth=0.7, zorder=4,
    )

    ax.set_xticks(x)
    ax.set_xticklabels(
        [pretty_solver_name(s) for s in solvers] if EMBED_AXIS_LABELS else [""] * len(solvers),
        rotation=30, ha="right",
    )

    if ylog:
        ax.set_yscale("log")

    _apply_labels(ax, ylabel=ylabel)
    ax.minorticks_on()

    savefig_all(fig, outpath)
    plt.close(fig)

def plot_box_by_solver(
        df: pd.DataFrame,
        solvers: list[str],
        colors: dict[str, str],
        ycol: str,
        ylabel: str,
        ylog: bool,
        figsize,
        outpath: Path,
):
    """
    Box plot per solver (distribution across instances).
    Median shown as a colored line; mean as a small diamond marker.
    """
    fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)

    data = []
    valid_solvers = []
    for s in solvers:
        vals = df.loc[df["solver"] == s, ycol].dropna().to_numpy(dtype=float)
        if len(vals) > 0:  # Only include solvers with data
            if ylog:
                vals = np.maximum(vals, LOG_EPS)
            data.append(vals)
            valid_solvers.append(s)

    # Skip if no valid data
    if len(data) == 0:
        print(f"DEBUG: No data for {ycol} - df shape: {df.shape}, ycol in df: {ycol in df.columns}, solvers: {solvers}")
        plt.close(fig)
        return

    bp = ax.boxplot(
        data,
        patch_artist=True,
        showfliers=False,
        widths=0.55,
        medianprops=dict(linewidth=1.6, color="black"),
        whiskerprops=dict(linewidth=0.7, linestyle="--", color="#555555"),
        capprops=dict(linewidth=0.7, color="#555555"),
        boxprops=dict(linewidth=0.7),
    )

    # Handle color assignment safely
    if bp.get("boxes"):
        for i, box in enumerate(bp["boxes"]):
            if i < len(valid_solvers):
                box.set_facecolor(colors[valid_solvers[i]])
                box.set_alpha(0.40)

    x = np.arange(1, len(valid_solvers) + 1)
    ax.set_xticks(x)
    ax.set_xticklabels(
        [pretty_solver_name(s) for s in valid_solvers] if EMBED_AXIS_LABELS else [""] * len(valid_solvers),
        rotation=30, ha="right",
    )

    if ylog:
        ax.set_yscale("log")
    _apply_labels(ax, ylabel=ylabel)
    ax.minorticks_on()


    savefig_all(fig, outpath)
    plt.close(fig)


def _fit_powerlaw_line(x: np.ndarray, y: np.ndarray):
    """
    Fit y ≈ a * x^b using log-log linear regression.
    Returns (a, b) or (None, None) if insufficient data.
    """
    mask = np.isfinite(x) & np.isfinite(y) & (x > 0) & (y > 0)
    x = x[mask]
    y = y[mask]
    if x.size < 3:
        return None, None
    lx = np.log10(x)
    ly = np.log10(y)
    b, loga = np.polyfit(lx, ly, deg=1)
    a = 10 ** loga
    return a, b


def plot_scatter_cloud(
        df: pd.DataFrame,
        solvers: list[str],
        colors: dict[str, str],
        markers: dict[str, str],
        xcol: str,
        ycol: str,
        xlabel: str,
        ylabel: str,
        xlog: bool,
        ylog: bool,
        figsize,
        outpath: Path,
        add_scaling_line: bool = True,
        linestyles: dict[str, str] | None = None,
):
    """
    Scatter plot where each dot is one instance.
    Optionally overlays a per-solver power-law scaling trend line.
    """
    fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)

    handles = []
    labels = []

    for s in solvers:
        sub = df[df["solver"] == s].copy()
        if sub.empty:
            continue

        x = pd.to_numeric(sub[xcol], errors="coerce").to_numpy(dtype=float)
        y = pd.to_numeric(sub[ycol], errors="coerce").to_numpy(dtype=float)

        if ylog:
            y = np.maximum(y, LOG_EPS)

        h = ax.scatter(
            x, y,
            marker=markers[s],
            s=5,
            alpha=0.55,
            color=colors[s],
            edgecolors="none",
            label=pretty_solver_name(s),
            zorder=3,
        )
        handles.append(h)
        labels.append(pretty_solver_name(s))

        if add_scaling_line and x.size >= 3:
            a, b = _fit_powerlaw_line(x, y)
            if a is not None:
                xpos = x[x > 0]
                x_min, x_max = np.nanmin(xpos), np.nanmax(xpos)
                if xlog:
                    xs = np.logspace(np.log10(x_min), np.log10(x_max), 200)
                else:
                    xs = np.linspace(x_min, x_max, 200)
                ys = a * (xs ** b)
                if ylog:
                    ys = np.maximum(ys, LOG_EPS)
                ls = linestyles[s] if linestyles else "--"
                ax.plot(xs, ys, linestyle=ls, linewidth=0.9, color=colors[s], alpha=0.85, zorder=2)

    if xlog:
        ax.set_xscale("log")
    if ylog:
        ax.set_yscale("log")

    _apply_labels(ax, xlabel, ylabel)
    ax.minorticks_on()

    add_legend(ax, handles, labels)

    savefig_all(fig, outpath)
    plt.close(fig)


def plot_runtime_decomposition_grouped_csv(
        df: pd.DataFrame,
        solvers: list[str],
        figsize,
        outpath: Path,
):
    required = ["transformation_time_micro_seconds", "solve_time_micro_seconds"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns for decomposition plot: {missing}")

    rows = []
    for s in solvers:
        sub = df[df["solver"] == s]
        if sub.empty:
            continue
        rows.append({
            "solver": s,
            "transform": float(np.median(sub["transformation_time_micro_seconds"].fillna(0))),
            "solve": float(np.median(sub["solve_time_micro_seconds"].fillna(0))),
        })

    if not rows:
        return

    agg = pd.DataFrame(rows).set_index("solver").reindex(solvers).dropna(how="all").reset_index()

    transform_vals = agg["transform"].fillna(0).to_numpy(dtype=float)
    solve_vals     = agg["solve"].fillna(0).to_numpy(dtype=float)
    totals         = transform_vals + solve_vals

    # Avoid division by zero for solvers where both are 0 / NaN
    safe_totals = np.where(totals > 0, totals, 1.0)
    pct_transform = transform_vals / safe_totals * 100.0
    pct_solve     = solve_vals     / safe_totals * 100.0

    fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)

    x     = np.arange(len(agg))
    width = 0.50

    # Fixed component colors — colorblind-safe (Wong palette)
    ax.bar(x, pct_transform, width,
           label="Transformation time",
           color="#D55E00", edgecolor="white", linewidth=0.4, alpha=0.88, zorder=3)
    ax.bar(x, pct_solve, width,
           bottom=pct_transform,
           label="Computation time",
           color="#009E73", edgecolor="white", linewidth=0.4, alpha=0.88, zorder=3)

    ax.set_ylim(0, 100)
    ax.yaxis.set_major_formatter(mticker.PercentFormatter())

    ax.set_xticks(x)
    ax.set_xticklabels(
        [pretty_solver_name(s) for s in agg["solver"]] if EMBED_AXIS_LABELS else [""] * len(agg),
        rotation=30, ha="right",
    )

    _apply_labels(ax, ylabel="Share of runtime [%]")

    savefig_all(fig, outpath)
    plt.close(fig)

def _graph_short_name(name: str) -> str:
    """Return a compact display name for a graph path/identifier.

    Always strips the file extension (e.g. .lgf, .gr, .dimacs) and any
    trailing run-index suffixes (_0, _1, …).
    """
    # Always take the last path component and drop the extension
    base = Path(name).stem
    # Remove trailing _N suffixes like _0, _1 that come from repeated runs
    base = re.sub(r"_\d+$", "", base)
    return base


def plot_metric_table_by_graph(
        df: pd.DataFrame,
        solvers: list[str],
        metric_col: str,
        metric_label: str,
        fmt: "Callable[[float], str]",
        outpath: Path,
        dedup_by_graph: bool = False,
):
    """
    Render a matplotlib table for a single metric where:
      • Rows    = solvers  (pretty names in the first column)
      • Columns = individual graph names  (sorted)
      • Cells   = median of *metric_col* for that (solver, graph) pair,
                  or "—" when no data is available.

    Parameters
    ----------
    dedup_by_graph : if True, drop duplicate (graph, solver) pairs before
                     aggregating (useful for oblivious_ratio which is unique
                     per graph rather than per run).
    """
    if metric_col not in df.columns or "graph" not in df.columns:
        return

    work = df[df["solver"].isin(solvers)][["solver", "graph", metric_col]].copy()
    work[metric_col] = pd.to_numeric(work[metric_col], errors="coerce")
    work = work.dropna(subset=[metric_col])

    if work.empty:
        return

    if dedup_by_graph:
        work = work.drop_duplicates(subset=["solver", "graph"])

    # Aggregate: median over repeated runs of the same (solver, graph)
    agg = (
        work.groupby(["solver", "graph"])[metric_col]
        .median()
        .reset_index()
    )

    graphs_unsorted = agg["graph"].unique()

    # ── Sort graphs by edge count (ascending) ─────────────────────────────────
    if "num_edges" in df.columns:
        edge_lookup = (
            df[["graph", "num_edges"]]
            .dropna(subset=["num_edges"])
            .groupby("graph")["num_edges"]
            .median()
        )
        graphs = sorted(
            graphs_unsorted,
            key=lambda g: (float(edge_lookup.get(g, np.inf)), _graph_short_name(g)),
        )
    else:
        graphs = sorted(graphs_unsorted, key=_graph_short_name)

    g_labels = [_graph_short_name(g) for g in graphs]

    # ── Build cell matrix (formatted strings + parallel raw floats) ──────────
    rows      = []   # formatted strings for display
    raw_vals  = []   # raw floats (np.nan for missing) — used to find minima
    active_solvers = []
    for s in solvers:
        sub = agg[agg["solver"] == s]
        if sub.empty:
            continue
        active_solvers.append(s)
        row     = [pretty_solver_name(s)]
        row_raw = [np.nan]   # solver name column — never a minimum
        for g in graphs:
            match = sub[sub["graph"] == g][metric_col]
            if not match.empty:
                v = float(match.iloc[0])
                row.append(fmt(v))
                row_raw.append(v)
            else:
                row.append("—")
                row_raw.append(np.nan)
        rows.append(row)
        raw_vals.append(row_raw)

    if not rows:
        return

    # ── Per-column minimum row index (data columns only, j >= 1) ─────────────
    n_rows = len(rows)
    n_cols = len(rows[0])   # == 1 (solver) + len(graphs)
    col_min_row: dict[int, int] = {}   # j -> row-index i of minimum
    for j in range(1, n_cols):
        col_vals = [raw_vals[i][j] for i in range(n_rows)]
        valid    = [(v, i) for i, v in enumerate(col_vals) if np.isfinite(v)]
        if valid:
            col_min_row[j] = min(valid, key=lambda t: t[0])[1]

    header  = ["Solver"] + g_labels

    # ── Figure sizing ─────────────────────────────────────────────────────────
    row_h      = 0.30
    solver_w   = 1.60   # first column (solver name) is wider
    data_col_w = 1.20   # data columns
    fig_h = row_h * (n_rows + 1) + 0.15
    fig_w = solver_w + data_col_w * (n_cols - 1)

    fig, ax = plt.subplots(figsize=(fig_w, fig_h), constrained_layout=True)
    ax.axis("off")

    tbl = ax.table(
        cellText=rows,
        colLabels=header,
        loc="center",
        cellLoc="center",
    )
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(7)

    # ── Header styling (B&W) ──────────────────────────────────────────────────
    for j in range(n_cols):
        cell = tbl[0, j]
        cell.set_facecolor("white")
        cell.set_text_props(color="black", fontweight="bold")
        cell.set_edgecolor("black")
        cell.set_linewidth(0.7)
        if j > 0:
            cell.get_text().set_rotation(45)
            cell.get_text().set_horizontalalignment("left")

    # ── Data-row styling (B&W) ────────────────────────────────────────────────
    for i in range(n_rows):
        bg = "#F2F2F2" if i % 2 == 0 else "white"
        for j in range(n_cols):
            cell = tbl[i + 1, j]
            cell.set_facecolor(bg)
            cell.set_edgecolor("#888888")
            cell.set_linewidth(0.4)
            if j == 0:
                cell.get_text().set_horizontalalignment("left")
            # Bold + slightly larger font for the minimum value in each column
            if j >= 1 and col_min_row.get(j) == i:
                cell.get_text().set_fontweight("bold")
                cell.get_text().set_fontsize(8)

    tbl.auto_set_column_width(list(range(n_cols)))

    savefig_all(fig, outpath)
    plt.close(fig)


def plot_error_lines_vs_edges(
        df: pd.DataFrame,
        solvers: list[str],
        colors: dict[str, str],
        figsize,
        outpath: Path,
        title: str | None = None,
        linestyles: dict[str, str] | None = None,
):
    """
    Line plot of mean relative error vs. number of edges, one line per solver.
    No individual scatter points — only the aggregated mean line and an
    optional shaded ±1 std band (controlled by SHOW_ERROR_BARS).
    """
    fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)

    handles = []
    labels = []

    for s in solvers:
        sub = df[df["solver"] == s][["num_edges", "relative_error"]].dropna()
        if sub.empty:
            continue

        agg = (
            sub.groupby("num_edges")["relative_error"]
            .agg(mean="mean", std="std")
            .reset_index()
            .sort_values("num_edges")
        )
        agg["std"] = agg["std"].fillna(0.0)

        x = agg["num_edges"].to_numpy(dtype=float)
        y = agg["mean"].to_numpy(dtype=float)
        e = agg["std"].to_numpy(dtype=float)

        ls = linestyles[s] if linestyles else "-"
        (h,) = ax.plot(x, y, linewidth=0.7, linestyle=ls, color=colors[s],
                       label=pretty_solver_name(s))
        handles.append(h)
        labels.append(pretty_solver_name(s))

        if SHOW_ERROR_BARS and np.any(e > 0):
            ax.fill_between(x, np.maximum(y - e, 0), y + e,
                            alpha=0.10, color=colors[s], linewidth=0)

    _apply_labels(ax, xlabel="Number of edges", ylabel="Mean relative error [%]", title=title)
    ax.minorticks_on()
    add_legend(ax, handles, labels)

    savefig_all(fig, outpath)
    plt.close(fig)


# ======================
# STANDALONE LEGEND PLATE
# ======================

def plot_legend_plate(
        solvers: list[str],
        colors: dict[str, str],
        markers: dict[str, str],
        linestyles: dict[str, str],
        outpath: Path,
        ncol: int | None = None,
):
    """
    Emit a figure that contains *only* the legend — no axes, no data.

    Each row shows the solver's line style + marker + colour + human-readable
    name.  Drop ``legend_plate.pdf`` once into a LaTeX paper and reference it
    from every figure caption.

    Parameters
    ----------
    ncol : number of legend columns.  Defaults to min(3, n_solvers).
    """
    n = len(solvers)
    if n == 0:
        return
    ncol = ncol or min(3, n)

    # Build proxy artists — a Line2D shows both the line style and the marker
    handles = []
    labels  = []
    for s in solvers:
        ls = linestyles.get(s, "-")
        h = mlines.Line2D(
            [], [],
            color=colors[s],
            linestyle=ls,
            linewidth=1.2,
            marker=markers[s],
            markersize=5,
            markeredgewidth=0.5,
            markeredgecolor="white",
            label=pretty_solver_name(s),
        )
        handles.append(h)
        labels.append(pretty_solver_name(s))

    # Size the figure to comfortably fit the legend
    nrows     = int(np.ceil(n / ncol))
    row_h_in  = 0.22   # height per legend row
    col_w_in  = 1.8    # width per legend column
    width_in  = col_w_in * ncol + 0.15
    height_in = row_h_in * nrows + 0.15

    fig, ax = plt.subplots(figsize=(width_in, height_in))
    ax.set_axis_off()
    fig.patch.set_facecolor("none")   # transparent background

    legend2 = ax.legend(
        handles, labels,
        ncol=ncol,
        loc="center",
        frameon=False,
        fontsize=8,
        handlelength=2.2,
        handleheight=0.9,
        handletextpad=0.5,
        columnspacing=1.2,
        labelspacing=0.4,
        borderpad=0.0,
    )

    fig.tight_layout(pad=0.0)
    savefig_all(fig, outpath)
    plt.close(fig)


# ======================
# QUALITY VS RUNTIME PLOT
# ======================

def plot_quality_vs_runtime(
        df: pd.DataFrame,
        solvers: list[str],
        colors: dict[str, str],
        markers: dict[str, str],
        figsize,
        outpath: Path,
        linestyles: dict[str, str] | None = None,
):
    """
    Scatter plot correlating solution quality against runtime.

    X-axis (log): total running time [microseconds] — one point per (graph, solver) instance.
    Y-axis (linear): oblivious ratio (lower = better; optimal = 1.0).

    A dashed y = 1 reference marks the theoretical optimum.

    Companion panel: per-solver mean ± 1 std summary so the quality/speed
    trade-off across solvers is immediately visible.
    """
    group_cols = ["solver", "graph"] if "graph" in df.columns else ["solver"]

    # One row per (graph, solver): mean total time + mean oblivious ratio.
    agg = (
        df[df["solver"].isin(solvers)]
        .groupby(group_cols, as_index=False)
        .agg(
            total_time_micro_seconds=("total_time_micro_seconds", "mean"),
            oblivious_ratio=("oblivious_ratio", "mean"),
        )
    )

    # Per-solver summary: mean & std of both axes
    summary = (
        agg.groupby("solver", as_index=False)
        .agg(
            t_mean=("total_time_micro_seconds",   "mean"),
            t_std= ("total_time_micro_seconds",   "std"),
            q_mean=("oblivious_ratio", "mean"),
            q_std= ("oblivious_ratio", "std"),
            n=     ("total_time_micro_seconds",   "count"),
        )
    )
    summary["t_std"] = summary["t_std"].fillna(0.0)
    summary["q_std"] = summary["q_std"].fillna(0.0)

    # ── Panel 1: individual instances ────────────────────────────────────────
    fig1, ax1 = plt.subplots(figsize=figsize, constrained_layout=True)

    ax1.axhline(1.0, color="#888888", linewidth=0.7, linestyle="--",
                zorder=1, label="_nolegend_")

    handles1, labels1 = [], []
    for s in solvers:
        sub = agg[agg["solver"] == s].copy()
        if sub.empty:
            continue
        x = pd.to_numeric(sub["total_time_micro_seconds"],   errors="coerce").to_numpy(dtype=float)
        y = pd.to_numeric(sub["oblivious_ratio"], errors="coerce").to_numpy(dtype=float)
        mask = np.isfinite(x) & np.isfinite(y) & (x > 0)
        x, y = x[mask], y[mask]
        if x.size == 0:
            continue

        h = ax1.scatter(
            x, y,
            marker=markers[s],
            s=10,
            alpha=0.65,
            color=colors[s],
            edgecolors="none",
            label=pretty_solver_name(s),
            zorder=3,
        )
        handles1.append(h)
        labels1.append(pretty_solver_name(s))

    ax1.set_xscale("log")
    _apply_labels(ax1, xlabel="Total running time [microseconds]", ylabel="Oblivious ratio")
    ax1.minorticks_on()
    add_legend(ax1, handles1, labels1)
    savefig_all(fig1, outpath.parent / (outpath.name + "_instances"))
    plt.close(fig1)

    # ── Panel 2: per-solver mean ± std (Pareto / summary view) ──────────────
    fig2, ax2 = plt.subplots(figsize=figsize, constrained_layout=True)

    ax2.axhline(1.0, color="#888888", linewidth=0.7, linestyle="--",
                zorder=1, label="_nolegend_")

    handles2, labels2 = [], []
    for s in solvers:
        row = summary[summary["solver"] == s]
        if row.empty:
            continue
        t_m  = float(row["t_mean"].iloc[0])
        t_s  = float(row["t_std"].iloc[0])
        q_m  = float(row["q_mean"].iloc[0])
        q_s  = float(row["q_std"].iloc[0])
        if not (np.isfinite(t_m) and t_m > 0 and np.isfinite(q_m)):
            continue

        ls = linestyles[s] if linestyles else "-"

        # cross-hairs: horizontal error bar (time std) + vertical (quality std)
        ax2.errorbar(
            t_m, q_m,
            xerr=[[min(t_s, t_m * 0.999)], [t_s]],   # log-safe lower bound
            yerr=[[min(q_s, max(q_m - 1e-9, 0))], [q_s]],
            fmt="none",
            ecolor=colors[s],
            elinewidth=0.7,
            capsize=2.5,
            capthick=0.7,
            alpha=0.55,
            zorder=2,
        )
        h = ax2.scatter(
            t_m, q_m,
            marker=markers[s],
            s=18,
            color=colors[s],
            edgecolors="white",
            linewidths=0.4,
            zorder=4,
            label=pretty_solver_name(s),
        )
        handles2.append(h)
        labels2.append(pretty_solver_name(s))

    ax2.set_xscale("log")
    _apply_labels(ax2, xlabel="Total running time [microseconds]", ylabel="Oblivious ratio (mean ± std)")
    ax2.minorticks_on()
    add_legend(ax2, handles2, labels2)
    savefig_all(fig2, outpath.parent / (outpath.name + "_summary"))
    plt.close(fig2)


# ======================
# MAIN PLOTTING ROUTINE (shared helper)
# ======================

def run_plots(df: pd.DataFrame, solvers: list, colors: dict, markers: dict,
              linestyles: dict, out_dir: Path, ylog: bool = False,
              mendel_mode: bool = False):
    """
    Generate the full set of plots for the given solver list into out_dir.

    Parameters
    ----------
    ylog        : if True all scaling/time axes use log-scale (used for Mendel pass).
    mendel_mode : if True also produce the Mendel-total-time plots.
    """
    out_dir.mkdir(parents=True, exist_ok=True)

    # ── Standalone legend plate (one file used by all figures in the paper) ──
    plot_legend_plate(
        solvers, colors, markers, linestyles,
        outpath=out_dir / "legend_plate",
    )

    _LP_SOLVERS = {"cohen", "lp"}
    solvers_mwu = [s for s in solvers if s not in _LP_SOLVERS]

    # ── Runtime scaling ──────────────────────────────────────────────────────
    agg_runtime = aggregate_mean_std(df[df["solver"].isin(solvers)].copy(), "total_time_micro_seconds")
    plot_lines(
        agg_runtime, solvers, colors, markers,
        xcol="num_edges", y_mean_col="mean", y_std_col="std",
        xlabel="Number of edges",
        ylabel="Total running time [microseconds]",
        xlog=True, ylog=True,
        figsize=FIGSIZE_SINGLE,
        outpath=out_dir / "runtime_lines_vs_edges",
        linestyles=linestyles,
    )

    # ── Solve time scaling (MWU solve time, excluding transformation) ─────────
    agg_solve = aggregate_mean_std(df[df["solver"].isin(solvers_mwu)].copy(), "solve_time_micro_seconds")
    plot_lines(
        agg_solve, solvers_mwu, colors, markers,
        xcol="num_edges", y_mean_col="mean", y_std_col="std",
        xlabel="Number of edges",
        ylabel="Solve time [microseconds]",
        xlog=True, ylog=True,
        figsize=FIGSIZE_SINGLE,
        outpath=out_dir / "solve_time_lines_vs_edges",
        linestyles=linestyles,
    )

    # ── MWU iterations ───────────────────────────────────────────────────────
    agg_mwu = aggregate_mean_std(df[df["solver"].isin(solvers_mwu)].copy(), "mwu_iterations")
    plot_lines(
        agg_mwu, solvers_mwu, colors, markers,
        xcol="num_edges", y_mean_col="mean", y_std_col="std",
        xlabel="Number of edges",
        ylabel="MWU iterations",
        xlog=True, ylog=ylog,
        figsize=FIGSIZE_SINGLE,
        outpath=out_dir / "mwu_iterations_lines_vs_edges",
        linestyles=linestyles,
    )

    # ── Oracle time ──────────────────────────────────────────────────────────
    agg_oracle = aggregate_mean_std(df[df["solver"].isin(solvers_mwu)].copy(), ORACLE_TIME_COL)
    plot_lines(
        agg_oracle, solvers_mwu, colors, markers,
        xcol="num_edges", y_mean_col="mean", y_std_col="std",
        xlabel="Number of edges",
        ylabel="Average oracle running time [microseconds]",
        xlog=True, ylog=True,
        figsize=FIGSIZE_SINGLE,
        outpath=out_dir / "oracle_time_lines_vs_edges",
        linestyles=linestyles,
    )

    # ── Load computation time (only for electrical flow) ────────────────────
    if "load_computation_micro_seconds" in df.columns:
        df_load = df[df["solver"].isin(solvers_mwu)].copy()
        # Only plot if there are non-NaN load computation values
        if df_load["load_computation_micro_seconds"].notna().any():
            agg_load = aggregate_mean_std(df_load, "load_computation_micro_seconds")
            if not agg_load.empty:
                plot_lines(
                    agg_load, solvers_mwu, colors, markers,
                    xcol="num_edges", y_mean_col="mean", y_std_col="std",
                    xlabel="Number of edges",
                    ylabel="Load computation time [microseconds]",
                    xlog=True, ylog=True,
                    figsize=FIGSIZE_SINGLE,
                    outpath=out_dir / "load_computation_time_lines_vs_edges",
                    linestyles=linestyles,
                )

    # ── Relative error box plots — always linear ─────────────────────────────
    plot_box_by_solver(
        df[df["solver"].isin(solvers)], solvers, colors,
        ycol="relative_error",
        ylabel="Relative error [%]",
        ylog=False,
        figsize=FIGSIZE_SINGLE,
        outpath=out_dir / "relative_error_box_by_solver",
    )

    if "demand_model" in df.columns:
        demand_models = sorted(df["demand_model"].dropna().unique())
        for demand in demand_models:
            df_demand = df[(df["demand_model"] == demand) & df["solver"].isin(solvers)]
            solvers_demand = [s for s in solvers
                              if not df_demand[df_demand["solver"] == s].empty]
            if not solvers_demand:
                continue
            plot_box_by_solver(
                df_demand, solvers_demand, colors,
                ycol="relative_error",
                ylabel="Relative error [%]",
                ylog=False,
                figsize=FIGSIZE_SINGLE,
                outpath=out_dir / f"relative_error_box_by_solver_{demand}",
            )

    # ── Transformation time ──────────────────────────────────────────────────
    df_t = df[df["solver"].isin(solvers_mwu)].copy()
    df_t["percent_transform_time"] = (
        df_t["transformation_time_micro_seconds"] / df_t["total_time_micro_seconds"] * 100.0
    )

    plot_box_by_solver(
        df_t, solvers_mwu, colors,
        ycol="percent_transform_time",
        ylabel="Transformation time in %",
        ylog=ylog,
        figsize=FIGSIZE_SINGLE,
        outpath=out_dir / "transformation_time_box_by_solver",
    )

    df_t["percent_transform_time"] = df_t["percent_transform_time"].replace(
        [np.inf, -np.inf], np.nan
    )
    df_t = df_t.dropna(subset=["percent_transform_time"])

    agg_transform_share = aggregate_mean_std_by_solver(df_t, "percent_transform_time")
    plot_solver_average_bars(
        agg_transform_share, solvers_mwu, colors,
        y_mean_col="mean", y_std_col="std",
        ylabel="Percentage of transformation time [%]",
        ylog=ylog,
        figsize=FIGSIZE_SINGLE,
        outpath=out_dir / "avg_transformation_time_share_by_solver",
    )

    plot_runtime_decomposition_grouped_csv(
        df[df["solver"].isin(solvers_mwu)], solvers_mwu,
        figsize=FIGSIZE_SINGLE,
        outpath=out_dir / "runtime_decomposition_grouped_csv",
    )

    # ── Per-metric tables: rows = solvers, columns = graph names ─────────────
    _df_mwu = df[df["solver"].isin(solvers_mwu)]

    plot_metric_table_by_graph(
        _df_mwu, solvers_mwu,
        metric_col="solve_time_micro_seconds",
        metric_label="Solve time [microseconds]",
        fmt=lambda v: f"{v:,.1f}",
        outpath=out_dir / "table_solve_time_by_graph",
    )
    plot_metric_table_by_graph(
        _df_mwu, solvers_mwu,
        metric_col="mwu_iterations",
        metric_label="MWU iterations",
        fmt=lambda v: f"{int(round(v)):,}",
        outpath=out_dir / "table_mwu_iterations_by_graph",
    )
    plot_metric_table_by_graph(
        _df_mwu, solvers_mwu,
        metric_col="avg_oracle_time_micro_seconds",
        metric_label="Avg oracle time [microseconds]",
        fmt=lambda v: f"{v:,.3f}",
        outpath=out_dir / "table_avg_oracle_time_by_graph",
    )
    plot_metric_table_by_graph(
        _df_mwu, solvers_mwu,
        metric_col="oblivious_ratio",
        metric_label="Oblivious ratio",
        fmt=lambda v: f"{v:.4f}",
        outpath=out_dir / "table_oblivious_ratio_by_graph",
        dedup_by_graph=True,
    )

    # ── Mendel total time (only in Mendel comparison pass) ───────────────────
    if mendel_mode and "mendel_total_micro_seconds" in df.columns:
        df_mendel = df[df["solver"].isin(solvers_mwu)].copy()
        df_mendel["mendel_total_micro_seconds"] = pd.to_numeric(
            df_mendel["mendel_total_micro_seconds"], errors="coerce"
        )
        df_mendel = df_mendel.dropna(subset=["mendel_total_micro_seconds"])
        # Only keep solvers that actually report mendel time (non-zero)
        solvers_mendel_time = [
            s for s in solvers_mwu
            if df_mendel.loc[df_mendel["solver"] == s, "mendel_total_micro_seconds"].gt(0).any()
        ]
        if solvers_mendel_time:
            agg_mendel = aggregate_mean_std(
                df_mendel[df_mendel["solver"].isin(solvers_mendel_time)], "mendel_total_micro_seconds"
            )
            plot_lines(
                agg_mendel, solvers_mendel_time, colors, markers,
                xcol="num_edges", y_mean_col="mean", y_std_col="std",
                xlabel="Number of edges",
                ylabel="Mendel scaling total time [microseconds]",
                xlog=False, ylog=ylog,
                figsize=FIGSIZE_SINGLE,
                outpath=out_dir / "mendel_total_time_lines_vs_edges",
                linestyles=linestyles,
            )
            plot_box_by_solver(
                df_mendel[df_mendel["solver"].isin(solvers_mendel_time)],
                solvers_mendel_time, colors,
                ycol="mendel_total_micro_seconds",
                ylabel="Mendel scaling total time [microseconds]",
                ylog=ylog,
                figsize=FIGSIZE_SINGLE,
                outpath=out_dir / "mendel_total_time_box_by_solver",
            )

    # ── Oblivious ratio ──────────────────────────────────────────────────────
    if "oblivious_ratio" in df.columns:
        df_oblivious = df[df["solver"].isin(solvers_mwu)].dropna(
            subset=["oblivious_ratio"]
        ).copy()
        df_oblivious["oblivious_ratio"] = pd.to_numeric(
            df_oblivious["oblivious_ratio"], errors="coerce"
        )
        df_oblivious = df_oblivious.dropna(subset=["oblivious_ratio"])

        dedup_cols = [c for c in ["graph", "solver", "num_edges", "oblivious_ratio"]
                      if c in df_oblivious.columns]
        df_oblivious_dedup = df_oblivious[dedup_cols].drop_duplicates(
            subset=["graph", "solver"] if "graph" in dedup_cols else None
        )

        solvers_oblivious = [s for s in solvers_mwu
                             if not df_oblivious_dedup[
                                 df_oblivious_dedup["solver"] == s].empty]

        if solvers_oblivious:
            plot_scatter_cloud(
                df_oblivious_dedup, solvers_oblivious, colors, markers,
                xcol="num_edges", ycol="oblivious_ratio",
                xlabel="Number of edges",
                ylabel="Oblivious ratio",
                xlog=False, ylog=False,
                figsize=FIGSIZE_SINGLE,
                outpath=out_dir / "oblivious_ratio_scatter_vs_edges",
                add_scaling_line=False,
                linestyles=linestyles,
            )

            plot_box_by_solver(
                df_oblivious_dedup, solvers_oblivious, colors,
                ycol="oblivious_ratio",
                ylabel="Oblivious ratio",
                ylog=False,
                figsize=FIGSIZE_SINGLE,
                outpath=out_dir / "oblivious_ratio_box_by_solver",
            )

            # ── Quality vs Runtime correlation ───────────────────────────────
            # Merge oblivious_ratio back onto rows that also have total runtime.
            qr_cols = ["solver", "total_time_micro_seconds", "oblivious_ratio"]
            if "graph" in df_oblivious.columns:
                qr_cols = ["solver", "graph", "total_time_micro_seconds", "oblivious_ratio"]
            df_qr = df_oblivious.loc[
                df_oblivious["solver"].isin(solvers_oblivious), qr_cols
            ].dropna(subset=["total_time_micro_seconds", "oblivious_ratio"]).copy()

            if not df_qr.empty:
                plot_quality_vs_runtime(
                    df_qr, solvers_oblivious, colors, markers,
                    figsize=FIGSIZE_SINGLE,
                    outpath=out_dir / "quality_vs_runtime",
                    linestyles=linestyles,
                )


# ======================
# MAIN
# ======================

def main():
    set_paper_style()

    df = pd.read_csv(RESULTS_CSV)
    df = normalize_time_columns_to_microseconds(df)

    if "status" in df.columns:
        # Only filter by status if there are non-NaN values
        if df["status"].notna().any():
            df = df[df["status"] == "OK"].copy()

    required_cols = {
        "solver", "num_edges", "total_time_micro_seconds", "mwu_iterations",
        "transformation_time_micro_seconds", "achieved_congestion", "offline_opt",
        ORACLE_TIME_COL,
    }
    missing = required_cols - set(df.columns)
    if missing:
        raise ValueError(f"Missing columns in results.csv: {sorted(missing)}")

    for c in ["num_edges", "total_time_micro_seconds", "solve_time_micro_seconds", "transformation_time_micro_seconds",
              "mwu_iterations", "achieved_congestion", "offline_opt", ORACLE_TIME_COL]:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")

    # CSV stores directed edge count; divide by 2 to get undirected edge count
    # used on all plot axes.
    df["num_edges"] = df["num_edges"] / 2

    # ── Derived columns ───────────────────────────────────────────────────────
    # 1) If solve_time is missing/NaN but MWU iterations and average oracle time
    #    are available, estimate solve_time = mwu_iterations × avg_oracle_time.
    if "solve_time_micro_seconds" in df.columns:
        if (ORACLE_TIME_COL in df.columns) and ("mwu_iterations" in df.columns):
            mask_missing_solve = df["solve_time_micro_seconds"].isna() | (df["solve_time_micro_seconds"] <= 0)
            if mask_missing_solve.any():
                df.loc[mask_missing_solve, "solve_time_micro_seconds"] = (
                    df.loc[mask_missing_solve, "mwu_iterations"]
                    * df.loc[mask_missing_solve, ORACLE_TIME_COL]
                )
                print(f"  ↳ Imputed {mask_missing_solve.sum()} missing solve_time_micro_seconds values "
                      f"from mwu_iterations × {ORACLE_TIME_COL}.")
    else:
        # Column not present at all — create it if we can
        if (ORACLE_TIME_COL in df.columns) and ("mwu_iterations" in df.columns):
            df["solve_time_micro_seconds"] = df["mwu_iterations"] * df[ORACLE_TIME_COL]
            print(f"  ↳ Created solve_time_micro_seconds from mwu_iterations × {ORACLE_TIME_COL}.")

    # 2) If total time is missing/NaN, fall back to solve time + transformation time.
    if "total_time_micro_seconds" in df.columns:
        mask_missing_total = df["total_time_micro_seconds"].isna()
        if mask_missing_total.any():
            for dep in ("solve_time_micro_seconds", "transformation_time_micro_seconds"):
                if dep not in df.columns:
                    df[dep] = np.nan
            df.loc[mask_missing_total, "total_time_micro_seconds"] = (
                df.loc[mask_missing_total, "solve_time_micro_seconds"].fillna(0.0)
                + df.loc[mask_missing_total, "transformation_time_micro_seconds"].fillna(0.0)
            )
            print(f"  ↳ Imputed {mask_missing_total.sum()} missing total_time_micro_seconds values "
                  f"from solve_time_micro_seconds + transformation_time_micro_seconds.")
    else:
        # Column not present at all
        for dep in ("solve_time_micro_seconds", "transformation_time_micro_seconds"):
            if dep not in df.columns:
                df[dep] = np.nan
        df["total_time_micro_seconds"] = (
            df["solve_time_micro_seconds"].fillna(0.0)
            + df["transformation_time_micro_seconds"].fillna(0.0)
        )
        print("  ↳ Created total_time_micro_seconds from solve_time_micro_seconds + transformation_time_micro_seconds.")

    df = df.dropna(subset=["solver", "num_edges"]).copy()

    df["ratio_pct"] = (
        df["achieved_congestion"] / df["offline_opt"].replace(0, np.nan)
    )
    df["relative_error"] = (df["ratio_pct"] ).clip(lower=0.0)

    print(f"DEBUG: After creating relative_error - df shape: {df.shape}, solvers in df: {df['solver'].nunique()}")

    # All solvers present in data
    # Filter out pointer HST variants (but keep flat HST and Pointer HST as separate entries)
    all_solvers = sorted(
        [s for s in df["solver"].unique() if "pointer" not in s.lower() or "HST)" in s],
        key=solver_sort_key,
    )
    print(f"DEBUG: all_solvers: {all_solvers}")
    print(f"DEBUG: Number of unique solvers: {len(all_solvers)}")

    # Shared visual style — same color/marker/linestyle per solver in every plot
    # IMPORTANT: Each solver variant (e.g., "Electrical Flow (naive)" vs "Electrical Flow (sketching)")
    # is a unique string and will get a unique color from the palette
    colors     = solver_palette_unique(all_solvers)
    markers    = solver_markers(all_solvers)
    linestyles = solver_linestyles(all_solvers)

    # Print color assignments for verification
    print("DEBUG: Color assignments:")
    for solver in all_solvers:
        print(f"  {solver:50s} → {colors[solver]}")

    # ── Solver sets ───────────────────────────────────────────────────────────
    # Pass 1: all solvers except Mendel/MendelScaling, naive electrical, and Pointer HST variants
    solvers_no_mendel = [s for s in all_solvers
                         if "Mendel" not in s and "MendelScaling" not in s
                         and "Electrical Flow (naive)" not in s
                         and "(Pointer HST)" not in s]

    # Pass 2: for each Mendel variant include its base + the mendel version
    # so the plot directly compares "before vs after" Mendel scaling.
    solvers_mendel_cmp: list[str] = []
    for s in all_solvers:
        if "Mendel" in s or "MendelScaling" in s:
            # Extract base name by removing the Mendel part
            base = s.replace(" + MendelScaling", "").replace(" (Mendel)", "")
            if base in all_solvers and base not in solvers_mendel_cmp:
                solvers_mendel_cmp.append(base)
            if s not in solvers_mendel_cmp:
                solvers_mendel_cmp.append(s)
    solvers_mendel_cmp = sorted(solvers_mendel_cmp, key=solver_sort_key)

    # Pass 3: for electrical flow variants, include both naive and sketching versions
    # to compare the two implementations
    solvers_electrical_cmp: list[str] = []
    for s in all_solvers:
        if "Electrical Flow" in s:
            solvers_electrical_cmp.append(s)
    solvers_electrical_cmp = sorted(solvers_electrical_cmp, key=solver_sort_key)

    # Pass 4: for HST data structures, include both Flat HST and Pointer HST versions
    # to compare the two tree implementations
    solvers_hst_cmp: list[str] = []
    for s in all_solvers:
        if "Raecke" in s or "Random MST" in s or "MST" in s:
            # Match Flat HST and Pointer HST pairs
            if "(Flat HST)" in s or "(Pointer HST)" in s:
                solvers_hst_cmp.append(s)
    solvers_hst_cmp = sorted(solvers_hst_cmp, key=solver_sort_key)

    # Pass 5: cycle removal strategy comparison, include only Flat HST tree-based solvers
    # (to compare naive vs Tarjan_SCC strategies)
    solvers_cycle_cmp: list[str] = []
    for s in all_solvers:
        if ("Raecke" in s or "Random MST" in s or "MST" in s) and "(Flat HST)" in s:
            solvers_cycle_cmp.append(s)
    solvers_cycle_cmp = sorted(solvers_cycle_cmp, key=solver_sort_key)

    # ── Pass 1: standard plots WITHOUT Mendel solvers (linear y-axes) ─────────
    print(f"[1/5] Generating no-Mendel plots (Flat HST & sketching variant only) → {OUT_DIR}")
    run_plots(df, solvers_no_mendel, colors, markers, linestyles, OUT_DIR,
              ylog=False, mendel_mode=False)

    # ── Pass 2: Mendel-comparison plots (log y-axes + Mendel time plots) ─────
    if solvers_mendel_cmp:
        mendel_dir = OUT_DIR / "mendel"
        print(f"[2/5] Generating Mendel-comparison plots → {mendel_dir}")
        run_plots(df, solvers_mendel_cmp, colors, markers, linestyles, mendel_dir,
                  ylog=True, mendel_mode=True)
    else:
        print("[2/5] No Mendel solvers found in data — skipping Mendel comparison plots.")

    # ── Pass 3: Electrical Flow variant comparison (naive vs sketching) ──────
    if solvers_electrical_cmp and len(solvers_electrical_cmp) > 1:
        electrical_dir = OUT_DIR / "electrical"
        print(f"[3/5] Generating Electrical Flow variant comparison plots → {electrical_dir}")
        run_plots(df, solvers_electrical_cmp, colors, markers, linestyles, electrical_dir,
                  ylog=False, mendel_mode=False)
    elif solvers_electrical_cmp:
        print("[3/5] Only one Electrical Flow variant found — skipping comparison plots.")
    else:
        print("[3/5] No Electrical Flow solvers found in data — skipping Electrical Flow comparison plots.")

    # ── Pass 4: HST data structure comparison (Flat HST vs Pointer HST) ──────
    if solvers_hst_cmp and len(solvers_hst_cmp) > len([s for s in solvers_hst_cmp if "(Flat HST)" in s]):
        hst_dir = OUT_DIR / "hst"
        print(f"[4/5] Generating HST data structure comparison plots → {hst_dir}")
        run_plots(df, solvers_hst_cmp, colors, markers, linestyles, hst_dir,
                  ylog=False, mendel_mode=False)
    else:
        print("[4/5] Only one HST data structure variant found — skipping HST comparison plots.")

    # ── Pass 5: Cycle removal strategy comparison (Flat HST only) ───────────
    if solvers_cycle_cmp:
        cycle_dir = OUT_DIR / "cycle_removal"
        print(f"[5/5] Generating cycle removal strategy comparison plots (Flat HST only) → {cycle_dir}")
        
        # Prepare data for cycle removal comparison: filter to Flat HST and naive/Tarjan_SCC only
        df_cycle = df.copy()
        
        # Filter to only Flat HST solvers
        df_cycle = df_cycle[df_cycle["solver"].isin(solvers_cycle_cmp)].copy()
        
        if "cycle_removal_type" in df_cycle.columns:
            # Filter to only naive and Tarjan_SCC strategies
            df_cycle = df_cycle[df_cycle["cycle_removal_type"].isin(["naive", "Tarjan_SCC"])].copy()
            
            # Create augmented solver names that include the strategy
            df_cycle["solver_with_strategy"] = (
                df_cycle["solver"] + " (" + df_cycle["cycle_removal_type"] + ")"
            )
            
            # Get unique solvers with strategies (only Flat HST)
            solvers_cycle_with_strategy = sorted(
                df_cycle["solver_with_strategy"].unique(),
                key=lambda s: solver_sort_key(s.split(" (")[0])  # Sort by base solver name
            )
            
            # Create color/marker/linestyle mappings for augmented solver names
            # Map augmented names back to base solver names for consistent colors
            colors_cycle = {}
            markers_cycle = {}
            linestyles_cycle = {}
            for augmented_name in solvers_cycle_with_strategy:
                base_name = augmented_name.rsplit(" (", 1)[0]  # Remove the (strategy) part
                if base_name in colors:
                    colors_cycle[augmented_name] = colors[base_name]
                    markers_cycle[augmented_name] = markers[base_name]
                    linestyles_cycle[augmented_name] = linestyles[base_name]
                else:
                    # Fallback to first color if base not found
                    colors_cycle[augmented_name] = _PALETTE[0]
                    markers_cycle[augmented_name] = ["o", "s", "^"][hash(augmented_name) % 3]
                    linestyles_cycle[augmented_name] = _LINESTYLES[hash(augmented_name) % len(_LINESTYLES)]
            
            # Temporarily rename solver column for plotting
            df_cycle_plot = df_cycle.copy()
            df_cycle_plot["solver"] = df_cycle_plot["solver_with_strategy"]
            
            # Generate plots with augmented solver names and augmented color mappings
            run_plots(df_cycle_plot, solvers_cycle_with_strategy, colors_cycle, markers_cycle, linestyles_cycle, cycle_dir,
                      ylog=False, mendel_mode=False)
        else:
            # Fallback if no cycle_removal_type column
            run_plots(df_cycle, solvers_cycle_cmp, colors, markers, 
                      linestyles, cycle_dir, ylog=False, mendel_mode=False)
    else:
        print("[5/5] No Flat HST tree-based solvers found in data — skipping cycle removal strategy plots.")

    print(f"✔ Paper-ready plots written to {OUT_DIR.resolve()}")
    print("✔ Unique solver colors (one per solver variant, consistent across all runs)")
    print("✔ Vector PDFs with embedded fonts (camera-ready)")
    print()
    print("Generated plot directories:")
    print(f"  [1] Main plots (Flat HST & sketching only):  {OUT_DIR.resolve()}")
    if solvers_mendel_cmp:
        print(f"  [2] Mendel comparison plots:                {(OUT_DIR / 'mendel').resolve()}")
    if solvers_electrical_cmp and len(solvers_electrical_cmp) > 1:
        print(f"  [3] Electrical Flow variant comparison:     {(OUT_DIR / 'electrical').resolve()}")
    if solvers_hst_cmp and len(solvers_hst_cmp) > len([s for s in solvers_hst_cmp if "(Flat HST)" in s]):
        print(f"  [4] HST data structure comparison:          {(OUT_DIR / 'hst').resolve()}")
    if solvers_cycle_cmp:
        print(f"  [5] Cycle removal strategy comparison:      {(OUT_DIR / 'cycle_removal').resolve()}")


if __name__ == "__main__":
    main()

