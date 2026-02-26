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
RESULTS_CSV = Path("/Users/halilibrahim/Desktop/Thesis/ObliviousRouting/results/combine.csv")
OUT_DIR = Path("/Users/halilibrahim/Desktop/Thesis/ObliviousRouting/plots/merged/")
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
        "electrical": "Electrical Flow",
        "electrical_parallel": "Electrical Flow (parallel)",
        "raecke_frt": "Räcke-FRT",
        "raecke_ckr": "Räcke–CKR",
        "raecke_mst": "Räcke–MST",
        "cohen": "LP (Applegate–Cohen)",
        "lp": "LP (Applegate–Cohen)",
        "random_mst": "Räcke–MST",
        "mst": "Räcke–MST",
        "frt": "Räcke-FRT",
        "ckr": "Räcke-CKR",
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

def plot_box_by_solver(
        df: pd.DataFrame,
        solvers: list[str],
        colors: dict[str, tuple],
        ycol: str,
        ylabel: str,
        ylog: bool,
        figsize,
        outpath: Path,
):
    """
    Box plot per solver (distribution across instances).
    Shows mean as a small marker on top of each box.
    """
    fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)

    data = []
    means = []
    for s in solvers:
        vals = df.loc[df["solver"] == s, ycol].dropna().to_numpy(dtype=float)
        if ylog:
            vals = np.maximum(vals, LOG_EPS)
        data.append(vals)
        means.append(np.mean(vals) if len(vals) > 0 else np.nan)

    bp = ax.boxplot(
        data,
        patch_artist=True,
        showfliers=False,     # publication-friendly; avoids extreme dot clutter
        widths=0.6,
        medianprops=dict(linewidth=0.9),
        whiskerprops=dict(linewidth=0.8),
        capprops=dict(linewidth=0.8),
    )

    for i, box in enumerate(bp["boxes"]):
        s = solvers[i]
        box.set_facecolor(colors[s])
        box.set_alpha(0.35)
        box.set_linewidth(0.8)

    # mean markers
    x = np.arange(1, len(solvers) + 1)
    ax.scatter(x, means, marker="D", s=14, zorder=3, color="black", linewidths=0.0, label="Mean")

    ax.set_xticks(x)
    ax.set_xticklabels([pretty_solver_name(s) for s in solvers], rotation=25, ha="right")


    ax.set_ylabel(ylabel)
    ax.minorticks_on()
    ax.grid(True, which="major")
    ax.grid(True, which="minor", alpha=0.10)

    # small legend just for the mean marker
    ax.legend(frameon=False, loc="upper right")

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
        colors: dict[str, tuple],
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
):
    """
    Scatter plot where each dot is one instance (no connecting lines).
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
            s=12,
            alpha=0.65,          # “cloud of dots”
            color=colors[s],
            edgecolors="none",
            label=pretty_solver_name(s),
        )
        handles.append(h)
        labels.append(pretty_solver_name(s))

        if add_scaling_line and x.size >= 3:
            a, b = _fit_powerlaw_line(x, y)
            if a is not None:
                xs = np.logspace(np.log10(np.nanmin(x[x > 0])), np.log10(np.nanmax(x)), 100)
                ys = a * (xs ** b)
                if ylog:
                    ys = np.maximum(ys, LOG_EPS)
                ax.plot(xs, ys, linestyle="--", linewidth=0.9, color=colors[s], alpha=0.9)

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


def plot_runtime_decomposition_grouped_csv(
        df: pd.DataFrame,
        solvers: list[str],
        figsize,
        outpath: Path,
):
    required = ["transformation_time_ms", "solve_time_ms", "total_time_ms"]
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
            "transform": float(np.median(sub["transformation_time_ms"])),
            "solve": float(np.median(sub["solve_time_ms"])),
            "total": float(np.median(sub["total_time_ms"])),
        })

    agg = pd.DataFrame(rows).set_index("solver").reindex(solvers).reset_index()

    fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)

    x = np.arange(len(agg))
    width = 0.24

    transform = np.maximum(agg["transform"].to_numpy(dtype=float), LOG_EPS)
    solve = np.maximum(agg["solve"].to_numpy(dtype=float), LOG_EPS)
    total = np.maximum(agg["total"].to_numpy(dtype=float), LOG_EPS)

    # Fixed component colors (super clear; works well in print)
    ax.bar(x - width, transform, width, label="Transformation time", color="#DD8452")
    ax.bar(x,         solve,     width, label="Computation time",          color="#55A868")
    ax.bar(x + width, total,     width, label="Total time",          color="#4C72B0")

    ax.set_xticks(x)
    ax.set_xticklabels([pretty_solver_name(s) for s in agg["solver"]], rotation=25, ha="right")

    ax.set_ylabel("Median time [ms]")
    ax.set_yscale("log")  # runtime plots usually much clearer on log scale

    ax.minorticks_on()
    ax.grid(True, which="major")
    ax.grid(True, which="minor", alpha=0.10)

    ax.legend(frameon=False, loc="upper left")

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
        "transformation_time_ms",
        "achieved_congestion",
        "offline_opt_value",
        ORACLE_TIME_COL,
    }
    missing = required_cols - set(df.columns)
    if missing:
        raise ValueError(f"Missing columns in results.csv: {sorted(missing)}")

    # Enforce numeric columns (prevents “ms”, “edges.” issues from silently breaking plots)
    for c in ["num_edges", "total_time_ms", "solve_time_ms", "transformation_time_ms",
              "mwu_iterations", "achieved_congestion", "offline_opt_value", ORACLE_TIME_COL]:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")

    df = df.dropna(subset=["solver", "num_edges"]).copy()
    df["relative_error"] = (abs(df["offline_opt_value"]-df["achieved_congestion"])) / df["achieved_congestion"]*100.0

    solvers = sorted(df["solver"].unique(), key=solver_sort_key)
    colors = solver_palette_unique(solvers)   # <-- UNIQUE COLORS
    markers = solver_markers(solvers)         # <-- helps grayscale

    # --- Scatter “cloud” plots (one dot per instance) ---
    plot_scatter_cloud(
        df, solvers, colors, markers,
        xcol="num_edges", ycol="total_time_ms",
        xlabel="Number of edges",
        ylabel="Total running time [ms]",
        xlog=False, ylog=True,
        figsize=FIGSIZE_SINGLE,
        outpath=OUT_DIR / "runtime_cloud_vs_edges",
        add_scaling_line=True,
    )

    plot_scatter_cloud(
        df, solvers, colors, markers,
        xcol="num_edges", ycol="mwu_iterations",
        xlabel="Number of edges",
        ylabel="MWU iterations",
        xlog=True, ylog=False,
        figsize=FIGSIZE_SINGLE,
        outpath=OUT_DIR / "mwu_iterations_cloud_vs_edges",
        add_scaling_line=True,
    )

    plot_scatter_cloud(
        df, solvers, colors, markers,
        xcol="num_edges", ycol=ORACLE_TIME_COL,
        xlabel="Number of edges",
        ylabel="Average oracle running time [ms]",
        xlog=True, ylog=True,
        figsize=FIGSIZE_SINGLE,
        outpath=OUT_DIR / "oracle_time_cloud_vs_edges",
        add_scaling_line=True,
    )

    # Bar plots of average metrics per solver
    # --- Box plots per solver (distribution across instances) ---
    plot_box_by_solver(
        df, solvers, colors,
        ycol="relative_error",
        ylabel="Relative error [%]",
        ylog=False,
        figsize=FIGSIZE_SINGLE,
        outpath=OUT_DIR / "relative_error_box_by_solver",
    )

    df["percent_transform_time"] = (
            df["transformation_time_ms"] / df["total_time_ms"] * 100.0
    )

    plot_box_by_solver(
        df, solvers, colors,
        ycol="percent_transform_time",
        ylabel="Transformation time [% of total time]",
        ylog=True,   # usually nicer; switch to False if you prefer linear
        figsize=FIGSIZE_SINGLE,
        outpath=OUT_DIR / "transformation_time_box_by_solver",
    )

    # Bar plots showing the percentage of solve time vs transformation time
    df["percent_transform_time"] = (
        df["transformation_time_ms"] / df["total_time_ms"] * 100.0
    )
    df["percent_transform_time"] = df["percent_transform_time"].replace([np.inf, -np.inf], np.nan)
    df = df.dropna(subset=["percent_transform_time"])

    agg_transform_share = aggregate_mean_std_by_solver(df, "percent_transform_time")

    plot_solver_average_bars(
        agg_transform_share, solvers, colors,
        y_mean_col="mean", y_std_col="std",
        ylabel="Percentage of transformation time [%]",
        ylog=False,
        figsize=FIGSIZE_SINGLE,
        outpath=OUT_DIR / "avg_transformation_time_share_by_solver",
    )


    plot_runtime_decomposition_grouped_csv(
        df, solvers,
        figsize=FIGSIZE_SINGLE,
        outpath=OUT_DIR / "runtime_decomposition_grouped_csv",
    )

    print(f"✔ Paper-ready plots written to {OUT_DIR.resolve()}")
    print("✔ Unique solver colors (no repeats up to 20 solvers; HSV beyond)")
    print("✔ Vector PDFs with embedded fonts (camera-ready)")


if __name__ == "__main__":
    main()
