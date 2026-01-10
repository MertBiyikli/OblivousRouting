#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt


# ======================
# CONFIG
# ======================
RESULTS_CSV = Path("/Users/halilibrahim/Desktop/Thesis/ObliviousRouting/results/results_SNDLIB_gravity.csv")
OUT_DIR = Path("/Users/halilibrahim/Desktop/Thesis/ObliviousRouting/plots/SNDLib/Gravity/")
OUT_DIR.mkdir(parents=True, exist_ok=True)

ORACLE_TIME_COL = "avg_oracle_time_ms"

# Typical LaTeX paper sizes
FIGSIZE_SINGLE = (3.35, 2.35)   # single-column
FIGSIZE_DOUBLE = (6.9, 3.0)     # double-column

EXPORT_PNG = True
DPI_RASTER = 450

SHOW_ERROR_BARS = False
LOG_EPS = 1e-3  # ms, safe lower bound for log-plots

# Legend control
LEGEND_NCOL = 2
LEGEND_OUTSIDE = True  # put legend above plot to avoid occluding data


# ======================
# STYLE (publication-grade)
# ======================

def set_paper_style():
    mpl.rcParams.update({
        # Typography
        "font.family": "serif",
        "font.serif": ["Times New Roman", "Times", "STIXGeneral", "DejaVu Serif"],
        "mathtext.fontset": "stix",

        "font.size": 7,
        "axes.labelsize": 7,
        "axes.titlesize": 7,
        "legend.fontsize": 6,
        "xtick.labelsize": 6,
        "ytick.labelsize": 6,

        # Lines/markers (slightly smaller = cleaner)
        "lines.linewidth": 0.75,
        "lines.markersize": 2.2,

        # Axes
        "axes.linewidth": 0.8,
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
        "xtick.minor.width": 0.6,
        "ytick.minor.width": 0.6,
        "xtick.direction": "in",
        "ytick.direction": "in",

        # Grid (subtle)
        "axes.grid": True,
        "grid.alpha": 0.22,
        "grid.linewidth": 0.55,

        # Output
        "pdf.fonttype": 42,  # embedded TrueType
        "ps.fonttype": 42,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.02,
    })


def savefig_all(fig: plt.Figure, outbase: Path):
    fig.savefig(outbase.with_suffix(".pdf"))
    if EXPORT_PNG:
        fig.savefig(outbase.with_suffix(".png"), dpi=DPI_RASTER)


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
    # Keep this conservative; adjust to your naming if you want
    mapping = {
        "electrical": "Electrical",
        "electrical_parallel": "Electrical (par)",
        "raecke_frt": "Räcke–FRT",
        "raecke_ckr": "Räcke–CKR",
        "cohen": "LP (Applegate–Cohen)",
        "random_mst": "Räcke–MST",
        "mst": "Räcke–MST",
    }
    return mapping.get(s, s)


def solver_sort_key(s: str) -> tuple:
    # Put "baseline/fast" first, then tree-based, then LP-ish, then others
    order = {
        "electrical": 0,
        "electrical_parallel": 1,
        "raecke_frt": 2,
        "raecke_ckr": 3,
        "random_mst": 4,
        "mst": 4,
        "cohen": 5,
    }
    return (order.get(s, 99), s)


# ======================
# COLORS (unique per solver, consistent across plots)
# ======================

def solver_palette_unique(solvers: list[str]) -> dict[str, tuple]:
    """
    Guarantees unique colors for all solvers (up to many solvers).
    - For <= 20: uses tab20 (high contrast).
    - For > 20: falls back to evenly spaced HSV colors.
    """
    n = len(solvers)
    colors: list[tuple]

    if n <= 20:
        cmap = plt.get_cmap("tab20")
        colors = [cmap(i) for i in range(n)]
    else:
        # evenly spaced hues
        cmap = plt.get_cmap("hsv")
        colors = [cmap(i / n) for i in range(n)]

    return {s: colors[i] for i, s in enumerate(solvers)}


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
    if not labels:
        return

    if LEGEND_OUTSIDE:
        ax.legend(
            handles, labels,
            loc="lower left",
            bbox_to_anchor=(0.0, 1.02),
            ncol=min(LEGEND_NCOL, max(1, len(labels))),
            frameon=False,
            borderaxespad=0.0,
            columnspacing=0.8,
            handletextpad=0.5,
        )
    else:
        ax.legend(frameon=False, ncol=min(LEGEND_NCOL, max(1, len(labels))))


def plot_lines(
        df_agg: pd.DataFrame,
        solvers: list[str],
        colors: dict[str, tuple],
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

        # Keep log-plots sane
        if ylog:
            y = np.maximum(y, LOG_EPS)
            e = np.minimum(e, y - LOG_EPS)

        h, = ax.plot(
            x, y,
            marker=markers[s],
            label=pretty_solver_name(s),
            color=colors[s],
        )
        handles.append(h)
        labels.append(pretty_solver_name(s))

        # optional uncertainty band (kept off by default)
        if SHOW_ERROR_BARS and np.any(e > 0):
            ax.fill_between(x, y - e, y + e, alpha=0.12, color=colors[s])

    if xlog:
        ax.set_xscale("log")
    if ylog:
        ax.set_yscale("log")

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    ax.minorticks_on()
    ax.grid(True, which="major")
    ax.grid(True, which="minor", alpha=0.10)

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
        colors: dict[str, tuple],
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
    bars = ax.bar(
        x, y,
        color=[colors[s] for s in solvers],
        alpha=0.8,
    )

    if SHOW_ERROR_BARS and np.any(e > 0):
        ax.errorbar(x, y, yerr=e, fmt="none", ecolor="black", capsize=2, linewidth=0.8)

    ax.set_xticks(x)
    ax.set_xticklabels([pretty_solver_name(s) for s in solvers], rotation=25, ha="right")

    if ylog:
        ax.set_yscale("log")

    ax.set_ylabel(ylabel)

    ax.minorticks_on()
    ax.grid(True, which="major")
    ax.grid(True, which="minor", alpha=0.10)

    savefig_all(fig, outpath)
    plt.close(fig)


# ======================
# MAIN
# ======================

def main():
    set_paper_style()

    df = pd.read_csv(RESULTS_CSV)

    # If you have status column from timeout script, keep only OK runs by default
    if "status" in df.columns:
        df = df[df["status"] == "OK"].copy()

    required_cols = {
        "solver",
        "num_edges",
        "total_time_ms",
        "mwu_iterations",
        "achieved_congestion",
        "offline_opt_value",
        ORACLE_TIME_COL,
    }
    missing = required_cols - set(df.columns)
    if missing:
        raise ValueError(f"Missing columns in results.csv: {sorted(missing)}")

    # Enforce numeric columns (prevents “ms”, “edges.” issues from silently breaking plots)
    for c in ["num_edges", "total_time_ms", "mwu_iterations", "achieved_congestion", "offline_opt_value", ORACLE_TIME_COL]:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    df = df.dropna(subset=["solver", "num_edges"]).copy()
    df["relative_error"] = (abs(df["offline_opt_value"]-df["achieved_congestion"])) / df["achieved_congestion"]*100.0

    solvers = sorted(df["solver"].unique(), key=solver_sort_key)
    colors = solver_palette_unique(solvers)   # <-- UNIQUE COLORS
    markers = solver_markers(solvers)         # <-- helps grayscale

    agg_runtime = aggregate_mean_std(df, "total_time_ms")
    agg_iters   = aggregate_mean_std(df, "mwu_iterations")
    agg_error   = aggregate_mean_std(df, "relative_error")
    agg_oracle  = aggregate_mean_std(df, ORACLE_TIME_COL)

    plot_lines(
        agg_runtime, solvers, colors, markers,
        "num_edges", "mean", "std",
        "Number of edges",
        "Total running time [ms]",
        xlog=True, ylog=True,
        figsize=FIGSIZE_SINGLE,
        outpath=OUT_DIR / "runtime_vs_edges",
    )

    plot_lines(
        agg_iters, solvers, colors, markers,
        "num_edges", "mean", "std",
        "Number of edges",
        "MWU iterations",
        xlog=True, ylog=False,
        figsize=FIGSIZE_SINGLE,
        outpath=OUT_DIR / "mwu_iterations_vs_edges",
    )

    plot_lines(
        agg_error, solvers, colors, markers,
        "num_edges", "mean", "std",
        "Number of edges",
        "Relative error (solver / optimal)",
        xlog=True, ylog=False,
        figsize=FIGSIZE_SINGLE,
        outpath=OUT_DIR / "error_vs_edges",
    )

    plot_lines(
        agg_oracle, solvers, colors, markers,
        "num_edges", "mean", "std",
        "Number of edges",
        "Average oracle running time [ms]",
        xlog=True, ylog=True,
        figsize=FIGSIZE_SINGLE,
        outpath=OUT_DIR / "avg_oracle_runtime_vs_edges",
    )


    avg_runtime_solver = aggregate_mean_std_by_solver(df, "total_time_ms")
    avg_error_solver   = aggregate_mean_std_by_solver(df, "relative_error")
    avg_oracle_solver  = aggregate_mean_std_by_solver(df, ORACLE_TIME_COL)
    avg_iters_solver   = aggregate_mean_std_by_solver(df, "mwu_iterations")

    plot_solver_average_bars(
        avg_runtime_solver, solvers, colors,
        y_mean_col="mean", y_std_col="std",
        ylabel="Average total running time [ms]",
        ylog=True,
        figsize=FIGSIZE_SINGLE,
        outpath=OUT_DIR / "avg_total_runtime_by_solver",
    )

    plot_solver_average_bars(
        avg_error_solver, solvers, colors,
        y_mean_col="mean", y_std_col="std",
        ylabel="Average relative error [%]",
        ylog=False,
        figsize=FIGSIZE_SINGLE,
        outpath=OUT_DIR / "avg_relative_error_by_solver",
    )

    plot_solver_average_bars(
        avg_iters_solver, solvers, colors,
        y_mean_col="mean", y_std_col="std",
        ylabel="Average MWU iterations",
        ylog=False,
        figsize=FIGSIZE_SINGLE,
        outpath=OUT_DIR / "avg_mwu_iterations_by_solver",
    )


    print(f"✔ Paper-ready plots written to {OUT_DIR.resolve()}")
    print("✔ Unique solver colors (no repeats up to 20 solvers; HSV beyond)")
    print("✔ Vector PDFs with embedded fonts (camera-ready)")


if __name__ == "__main__":
    main()
