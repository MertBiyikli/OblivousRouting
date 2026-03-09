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
RESULTS_CSV = Path("/Users/halilibrahim/Desktop/Thesis/ObliviousRouting/new_results/SNDLib/SNDLib.csv")
OUT_DIR = Path("/Users/halilibrahim/Desktop/Thesis/ObliviousRouting/plots/SNDLib/")
OUT_DIR.mkdir(parents=True, exist_ok=True)

ORACLE_TIME_COL = "avg_oracle_time_ms"

# Typical LaTeX paper sizes
FIGSIZE_SINGLE = (2.87, 2.17)   # thesis single-column-ish
FIGSIZE_DOUBLE = (5.91, 2.17)     # thesis double-column-ish

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
        "lines.linewidth": 1.2,
        "lines.markersize": 4.5,
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
    """
    order = {
        "electrical":          0,
        "electrical_parallel": 1,
        "raecke_ckr":          2,
        "ckr":                 2,
        "raecke_ckr_mendel":   3,
        "ckr_mendel":          3,
        "raecke_frt":          4,
        "frt":                 4,
        "raecke_frt_mendel":   5,
        "frt_mendel":          5,
        "raecke_mst":          6,
        "raecke_random_mst":   6,
        "random_mst":          6,
        "mst":                 6,
        "raecke_mst_mendel":   6,
        # LP variants always last
        "cohen":               99,
        "lp":                  99,
    }
    return (order.get(s, 90), s)


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
]

_LINESTYLES = ["-", "--", "-.", ":", (0, (3, 1, 1, 1)), (0, (5, 1))]


def solver_palette_unique(solvers: list[str]) -> dict[str, str]:
    """Unique color per solver; cycles through colorblind-safe palette."""
    return {s: _PALETTE[i % len(_PALETTE)] for i, s in enumerate(solvers)}


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
    if not labels:
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
        # space markers so they don't clutter when many points
        every = max(1, len(x) // 6)
        h, = ax.plot(
            x, y,
            linestyle=ls,
            marker=markers[s],
            markevery=every,
            markersize=4.5,
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

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
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
    ax.set_xticklabels([pretty_solver_name(s) for s in solvers], rotation=30, ha="right")

    if ylog:
        ax.set_yscale("log")

    ax.set_ylabel(ylabel)
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
    means = []
    for s in solvers:
        vals = df.loc[df["solver"] == s, ycol].dropna().to_numpy(dtype=float)
        if ylog:
            vals = np.maximum(vals, LOG_EPS)
        data.append(vals)
        means.append(float(np.mean(vals)) if len(vals) > 0 else np.nan)

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

    for i, box in enumerate(bp["boxes"]):
        box.set_facecolor(colors[solvers[i]])
        box.set_alpha(0.40)

    # mean markers
    x = np.arange(1, len(solvers) + 1)
    ax.scatter(
        x, means,
        marker="D", s=18, zorder=4,
        color="#222222", linewidths=0.0,
        label="Mean",
    )

    ax.set_xticks(x)
    ax.set_xticklabels([pretty_solver_name(s) for s in solvers], rotation=30, ha="right")

    if ylog:
        ax.set_yscale("log")
    ax.set_ylabel(ylabel)
    ax.minorticks_on()

    add_legend(ax, [ax.collections[-1]], ["Mean"])

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

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
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

    # Fixed component colors — colorblind-safe (Wong palette)
    ax.bar(x - width, transform, width, label="Transformation time",
           color="#D55E00", edgecolor="white", linewidth=0.4, alpha=0.88, zorder=3)
    ax.bar(x,         solve,     width, label="Computation time",
           color="#009E73", edgecolor="white", linewidth=0.4, alpha=0.88, zorder=3)
    ax.bar(x + width, total,     width, label="Total time",
           color="#0072B2", edgecolor="white", linewidth=0.4, alpha=0.88, zorder=3)

    ax.set_xticks(x)
    ax.set_xticklabels([pretty_solver_name(s) for s in agg["solver"]], rotation=30, ha="right")

    ax.set_ylabel("Median time [ms]")
    ax.set_yscale("log")

    ax.minorticks_on()

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

    ax.set_xlabel("Number of edges")
    ax.set_ylabel("Mean relative error [%]")
    if title:
        ax.set_title(title)

    ax.minorticks_on()
    add_legend(ax, handles, labels)

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
        "offline_opt",
        ORACLE_TIME_COL,
    }
    missing = required_cols - set(df.columns)
    if missing:
        raise ValueError(f"Missing columns in results.csv: {sorted(missing)}")

    # Enforce numeric columns (prevents “ms”, “edges.” issues from silently breaking plots)
    for c in ["num_edges", "total_time_ms", "solve_time_ms", "transformation_time_ms",
              "mwu_iterations", "achieved_congestion", "offline_opt", ORACLE_TIME_COL]:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")

    df = df.dropna(subset=["solver", "num_edges"]).copy()

    # ratio_pct = oblivious ratio as a percentage: how many % of the offline optimum
    # the scheme achieves. Always >= 100% (scheme can never beat offline optimum).
    df["ratio_pct"] = (
        df["achieved_congestion"] / df["offline_opt"].replace(0, np.nan) * 100.0
    )

    # relative_error = excess over the offline optimum in percent.
    # = (achieved - offline) / offline * 100  = ratio_pct - 100
    # Clamp to 0 to absorb floating-point noise where achieved ≈ offline.
    df["relative_error"] = (df["ratio_pct"] - 100.0).clip(lower=0.0)

    solvers = sorted(df["solver"].unique(), key=solver_sort_key)
    colors = solver_palette_unique(solvers)
    markers = solver_markers(solvers)
    linestyles = solver_linestyles(solvers)

    # --- Scaling plots: aggregated mean lines (cleaner than scatter clouds) ---
    agg_runtime = aggregate_mean_std(df, "total_time_ms")
    plot_lines(
        agg_runtime, solvers, colors, markers,
        xcol="num_edges", y_mean_col="mean", y_std_col="std",
        xlabel="Number of edges",
        ylabel="Total running time [ms]",
        xlog=False, ylog=True,
        figsize=FIGSIZE_SINGLE,
        outpath=OUT_DIR / "runtime_lines_vs_edges",
        linestyles=linestyles,
    )

    # Exclude LP-type solvers: they have no MWU iterations or oracle time
    _LP_SOLVERS = {"cohen", "lp"}
    solvers_mwu = [s for s in solvers if s not in _LP_SOLVERS]

    agg_mwu = aggregate_mean_std(df[df["solver"].isin(solvers_mwu)].copy(), "mwu_iterations")
    plot_lines(
        agg_mwu, solvers_mwu, colors, markers,
        xcol="num_edges", y_mean_col="mean", y_std_col="std",
        xlabel="Number of edges",
        ylabel="MWU iterations",
        xlog=True, ylog=False,
        figsize=FIGSIZE_SINGLE,
        outpath=OUT_DIR / "mwu_iterations_lines_vs_edges",
        linestyles=linestyles,
    )

    agg_oracle = aggregate_mean_std(df[df["solver"].isin(solvers_mwu)].copy(), ORACLE_TIME_COL)
    plot_lines(
        agg_oracle, solvers_mwu, colors, markers,
        xcol="num_edges", y_mean_col="mean", y_std_col="std",
        xlabel="Number of edges",
        ylabel="Average oracle running time [ms]",
        xlog=False, ylog=True,
        figsize=FIGSIZE_SINGLE,
        outpath=OUT_DIR / "oracle_time_lines_vs_edges",
        linestyles=linestyles,
    )

    # Bar plots of average metrics per solver
    # Bar plots per solver (distribution across instances) ---
    plot_box_by_solver(
        df, solvers, colors,
        ycol="relative_error",
        ylabel="Relative error [%]",
        ylog=False,
        figsize=FIGSIZE_SINGLE,
        outpath=OUT_DIR / "relative_error_box_by_solver",
    )

    # Per-demand-model relative error box plots
    if "demand_model" in df.columns:
        demand_models = sorted(df["demand_model"].dropna().unique())
        for demand in demand_models:
            df_demand = df[df["demand_model"] == demand]
            # Only keep solvers that have data for this demand model
            solvers_demand = [s for s in solvers if not df_demand[df_demand["solver"] == s].empty]
            if not solvers_demand:
                continue
            plot_box_by_solver(
                df_demand, solvers_demand, colors,
                ycol="relative_error",
                ylabel="Relative error [%]",
                ylog=False,
                figsize=FIGSIZE_SINGLE,
                outpath=OUT_DIR / f"relative_error_box_by_solver_{demand}",
            )


    df["percent_transform_time"] = (
            df["transformation_time_ms"] / df["total_time_ms"] * 100.0
    )

    plot_box_by_solver(
        df, solvers_mwu, colors,
        ycol="percent_transform_time",
        ylabel="Transformation time in %",
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
        agg_transform_share, solvers_mwu, colors,
        y_mean_col="mean", y_std_col="std",
        ylabel="Percentage of transformation time [%]",
        ylog=False,
        figsize=FIGSIZE_SINGLE,
        outpath=OUT_DIR / "avg_transformation_time_share_by_solver",
    )


    plot_runtime_decomposition_grouped_csv(
        df, solvers_mwu,
        figsize=FIGSIZE_SINGLE,
        outpath=OUT_DIR / "runtime_decomposition_grouped_csv",
    )

    print(f"✔ Paper-ready plots written to {OUT_DIR.resolve()}")
    print("✔ Unique solver colors (no repeats up to 20 solvers; HSV beyond)")
    print("✔ Vector PDFs with embedded fonts (camera-ready)")


if __name__ == "__main__":
    main()