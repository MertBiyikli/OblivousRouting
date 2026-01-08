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

RESULTS_CSV = Path("results/results_rocketfuel.csv")
OUT_DIR = Path("plots")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Fixed column name (matches stdout + CSV schema)
ORACLE_TIME_COL = "avg_oracle_time_ms"

# Typical LaTeX paper sizes
FIGSIZE_SINGLE = (3.35, 2.35)   # single-column
FIGSIZE_DOUBLE = (6.9, 3.0)     # double-column

EXPORT_PNG = True
DPI_RASTER = 350
SHOW_ERROR_BARS = False
LOG_EPS = 1e-3  # ms, safe lower bound for oracle time


# ======================
# STYLE (publication-grade)
# ======================

def set_paper_style():
    mpl.rcParams.update({
        # ---------- Typography ----------
        "font.family": "serif",
        "font.serif": [
            "Times New Roman",
            "Times",
            "STIXGeneral",
            "DejaVu Serif",
        ],
        "mathtext.fontset": "stix",

        "font.size": 7,
        "axes.labelsize": 7,
        "axes.titlesize": 7,
        "legend.fontsize": 6,
        "xtick.labelsize": 6,
        "ytick.labelsize": 6,

        # ---------- Lines ----------
        "lines.linewidth": 0.7,
        "lines.markersize": 2.5,

        # ---------- Axes ----------
        "axes.linewidth": 0.8,
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
        "xtick.minor.width": 0.6,
        "ytick.minor.width": 0.6,
        "xtick.direction": "in",
        "ytick.direction": "in",

        # ---------- Grid ----------
        "axes.grid": True,
        "grid.alpha": 0.25,
        "grid.linewidth": 0.6,

        # ---------- PDF output ----------
        "pdf.fonttype": 42,   # embedded TrueType (editable)
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
    """
    Mean / std aggregation per (solver, num_edges)
    """
    g = (
        df.groupby(["solver", "num_edges"], as_index=False)[ycol]
        .agg(["mean", "std", "count"])
        .reset_index()
    )
    g.rename(columns={"count": "n"}, inplace=True)
    g["std"] = g["std"].fillna(0.0)
    return g


# ======================
# PLOTTING HELPERS
# ======================

def solver_palette(solvers: list[str]) -> dict[str, tuple]:
    cmap = plt.get_cmap("tab10")
    return {s: cmap(i % 10) for i, s in enumerate(solvers)}


def plot_with_band(
        df_agg: pd.DataFrame,
        solvers: list[str],
        colors: dict[str, tuple],
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

    for s in solvers:
        sub = df_agg[df_agg["solver"] == s].sort_values(xcol)
        x = sub[xcol].to_numpy()
        y = sub[y_mean_col].to_numpy()
        e = sub[y_std_col].to_numpy()
        y = np.maximum(y, LOG_EPS)
        e = np.minimum(e, y - LOG_EPS)   # prevent y - e <= 0


        ax.plot(x, y, marker="o", label=s, color=colors[s])
        if SHOW_ERROR_BARS:
            ax.fill_between(x, y - e, y + e, alpha=0.15, color=colors[s])

    if xlog:
        ax.set_xscale("log")
    if ylog:
        ax.set_yscale("log")

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    ax.minorticks_on()
    ax.grid(True, which="major")
    ax.grid(True, which="minor", alpha=0.12)

    ax.legend(frameon=False)
    savefig_all(fig, outpath)
    plt.close(fig)


# ======================
# MAIN
# ======================

def main():
    set_paper_style()

    df = pd.read_csv(RESULTS_CSV)

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

    df = df.copy()
    df["relative_error"] = df["achieved_congestion"] / df["offline_opt_value"]

    solvers = sorted(df["solver"].unique())
    colors = solver_palette(solvers)

    agg_runtime = aggregate_mean_std(df, "total_time_ms")
    agg_iters   = aggregate_mean_std(df, "mwu_iterations")
    agg_error   = aggregate_mean_std(df, "relative_error")
    agg_oracle  = aggregate_mean_std(df, ORACLE_TIME_COL)

    # 1) Total runtime
    plot_with_band(
        agg_runtime, solvers, colors,
        "num_edges", "mean", "std",
        "Number of edges",
        "Total running time [ms]",
        xlog=True, ylog=True,
        figsize=FIGSIZE_SINGLE,
        outpath=OUT_DIR / "runtime_vs_edges",
    )

    # 2) MWU iterations
    plot_with_band(
        agg_iters, solvers, colors,
        "num_edges", "mean", "std",
        "Number of edges",
        "MWU iterations",
        xlog=True, ylog=False,
        figsize=FIGSIZE_SINGLE,
        outpath=OUT_DIR / "mwu_iterations_vs_edges",
    )

    # 3) Relative error
    plot_with_band(
        agg_error, solvers, colors,
        "num_edges", "mean", "std",
        "Number of edges",
        "Relative error (solver / optimal)",
        xlog=True, ylog=False,
        figsize=FIGSIZE_SINGLE,
        outpath=OUT_DIR / "error_vs_edges",
    )

    # 4) Average oracle runtime
    plot_with_band(
        agg_oracle, solvers, colors,
        "num_edges", "mean", "std",
        "Number of edges",
        "Average oracle running time [ms]",
        xlog=True, ylog=True,
        figsize=FIGSIZE_SINGLE,
        outpath=OUT_DIR / "avg_oracle_runtime_vs_edges",
    )

    print(f"✔ Paper-ready plots written to {OUT_DIR.resolve()}")
    print("✔ Vector PDFs with embedded fonts (camera-ready)")


if __name__ == "__main__":
    main()
