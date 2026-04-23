import re
import argparse
from collections.abc import Callable
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker


LOG_EPS = 1e-3  # microseconds, safe lower bound for log-plots
FIGSIZE_SINGLE = (2.87, 2.17)   # thesis single-column-ish
FIGSIZE_DOUBLE = (5.91, 2.17)     # thesis double-column-ish
def savefig_all(fig: plt.Figure, outbase: Path):
    """Save figure to PDF and optionally PNG with optimizations for speed."""
    # Save PDF (fastest vector format)
    fig.savefig(outbase.with_suffix(".pdf"), bbox_inches='tight', pad_inches=0.05)


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

    for s in solvers:
        sub = df_agg[df_agg["solver"] == s].sort_values(xcol)
        if sub.empty:
            continue

        x = pd.to_numeric(sub[xcol], errors="coerce").to_numpy(dtype=float)
        y = pd.to_numeric(sub[y_mean_col], errors="coerce").to_numpy(dtype=float)
        e = pd.to_numeric(sub[y_std_col], errors="coerce").to_numpy(dtype=float)

        mask = np.isfinite(x) & np.isfinite(y) & np.isfinite(e)
        if xlog:
            mask &= (x > 0)
        if ylog:
            mask &= (y > 0)

        x = x[mask]
        y = y[mask]
        e = e[mask]

        if len(x) == 0:
            continue

        if ylog:
            y = np.maximum(y, LOG_EPS)
            e = np.minimum(e, np.maximum(y - LOG_EPS, 0))

        # lighter confidence band instead of visually dominant error bars
        if len(e) > 0 and np.any(e > 0):
            y_low = np.maximum(y - e, LOG_EPS if ylog else -np.inf)
            y_high = y + e
            ax.fill_between(
                x, y_low, y_high,
                color=colors[s],
                alpha=0.08,
                linewidth=0,
                zorder=1,
            )

        # thin, slightly transparent line
        ax.plot(
            x, y,
            linestyle=linestyles[s] if linestyles else "-",
            linewidth=0.7,
            color=colors[s],
            alpha=0.7,
            zorder=2,
        )

        # emphasize points more than line
        ax.scatter(
            x, y,
            marker=markers[s],
            s=14,
            color=colors[s],
            edgecolors="white",
            linewidths=0.25,
            alpha=0.9,
            zorder=3,
        )

    if xlog:
        ax.set_xscale("log")
    if ylog:
        ax.set_yscale("log")

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.minorticks_on()

    if xlog:
        ax.xaxis.set_major_locator(mticker.LogLocator(base=10, numticks=6))
    else:
        ax.xaxis.set_major_locator(mticker.MaxNLocator(nbins=5))

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

    If ylog is False but values span more than 2 orders of magnitude,
    automatically enable logarithmic scaling for better visualization.

    IMPORTANT: Only includes instances where ALL solvers have data, to ensure
    a fair comparison without distribution shifts from missing data.
    """
    fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)

    # ── Filter to instances present in ALL solvers ───────────────────────────
    # Identify which instances (rows) have valid data for each solver
    if "graph" in df.columns:
        # Group by graph to identify instances that have data for all solvers
        instance_col = "graph"
    else:
        # If no graph column, use row index as instance identifier
        instance_col = None

    # Collect instances that have data for each solver
    valid_instances_per_solver = {}
    for s in solvers:
        df_s = df.loc[df["solver"] == s]
        # Get instances (graphs or indices) that have non-null values for ycol
        valid_rows = df_s[df_s[ycol].notna()]
        if instance_col and instance_col in valid_rows.columns:
            valid_instances_per_solver[s] = set(valid_rows[instance_col].unique())
        else:
            # Fall back to row indices
            valid_instances_per_solver[s] = set(valid_rows.index)

    # Find instances common to all solvers
    if valid_instances_per_solver:
        common_instances = set.intersection(*valid_instances_per_solver.values())
        # Filter dataframe to only rows with instances in common_instances
        if instance_col and instance_col in df.columns:
            df_filtered = df[df[instance_col].isin(common_instances)].copy()
        else:
            df_filtered = df.loc[list(common_instances)].copy()
    else:
        df_filtered = df.copy()

    data = []
    valid_solvers = []
    all_vals_for_range = []

    for s in solvers:
        vals = df_filtered.loc[df_filtered["solver"] == s, ycol].dropna().to_numpy(dtype=float)
        if len(vals) > 0:  # Only include solvers with data
            data.append(vals)
            valid_solvers.append(s)
            all_vals_for_range.extend(vals)

    # Skip if no valid data
    if len(data) == 0:
        print(f"DEBUG: No data for {ycol} - df shape: {df.shape}, ycol in df: {ycol in df.columns}, solvers: {solvers}")
        plt.close(fig)
        return

    # Auto-detect if log scale is needed
    actual_ylog = ylog
    if not ylog and len(all_vals_for_range) > 0:
        all_vals_array = np.array(all_vals_for_range)
        # Filter out zero and negative values for range calculation
        positive_vals = all_vals_array[all_vals_array > 0]
        if len(positive_vals) > 0:
            val_min = np.min(positive_vals)
            val_max = np.max(positive_vals)
            # If range spans more than 2 orders of magnitude (100x), use log scale
            if val_max / val_min > 100:
                actual_ylog = True

    # Apply log transformation if needed
    if actual_ylog:
        data = [np.maximum(d, LOG_EPS) for d in data]

    bp = ax.boxplot(
        data,
        patch_artist=True,
        showfliers=False,
        whis=[0, 100],
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

    if actual_ylog:
        ax.set_yscale("log")

    ax.minorticks_on()

    savefig_all(fig, outpath)
    plt.close(fig)

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
            markersize=2,
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

def export_electrical_sketching_summary_table(
        df: pd.DataFrame,
        out_dir: Path,
):
    """
    Create a compact comparison table for Electrical Flow (naive) vs
    Electrical Flow (sketching).

    Output is aggregated into 3 fixed edge-range buckets:
    [1, 500], [500, 15000], [15000, ...)

    Columns:
    - Number of edges
    - Naive [s]
    - Sketching [s]
    - Speedup
    - Delta oblivious ratio
    """
    naive_df = df[df["solver"] == "Electrical Flow (naive)"].copy()
    sketch_df = df[df["solver"] == "Electrical Flow (sketching)"].copy()

    if naive_df.empty or sketch_df.empty:
        print("WARNING: Missing naive or sketching electrical flow data")
        return

    # Ensure numeric columns
    for frame in (naive_df, sketch_df):
        frame["num_edges"] = pd.to_numeric(frame["num_edges"], errors="coerce")
        frame["total_time_micro_seconds"] = pd.to_numeric(
            frame["total_time_micro_seconds"], errors="coerce"
        )
        if "oblivious_ratio" in frame.columns:
            frame["oblivious_ratio"] = pd.to_numeric(
                frame["oblivious_ratio"], errors="coerce"
            )

    naive_grouped = (
        naive_df.groupby("num_edges").agg(
            naive_mean_us=("total_time_micro_seconds", "mean"),
            naive_obl_mean=("oblivious_ratio", "mean"),
        ).reset_index()
    )

    sketch_grouped = (
        sketch_df.groupby("num_edges").agg(
            sketch_mean_us=("total_time_micro_seconds", "mean"),
            sketch_obl_mean=("oblivious_ratio", "mean"),
        ).reset_index()
    )

    summary = naive_grouped.merge(sketch_grouped, on="num_edges", how="inner")

    if summary.empty:
        print("WARNING: No overlapping num_edges values between naive and sketching runs")
        return

    summary["naive_mean_s"] = summary["naive_mean_us"] / 1e6
    summary["sketch_mean_s"] = summary["sketch_mean_us"] / 1e6
    summary["speedup"] = summary["naive_mean_us"] / summary["sketch_mean_us"]
    summary["obl_ratio_diff"] = summary["sketch_obl_mean"] - summary["naive_obl_mean"]

    # Fixed buckets to match your desired table style
    bins = [1, 500, 15000, np.inf]
    labels = [r"$[1, 500]$", r"$[500,15000]$", r"$[15000, \cdots]$"]

    summary = summary[summary["num_edges"] >= 1].copy()
    summary["edge_bucket"] = pd.cut(
        summary["num_edges"],
        bins=bins,
        labels=labels,
        include_lowest=True,
        right=True,
    )

    compact = (
        summary.groupby("edge_bucket", observed=True)
        .agg(
            naive_s=("naive_mean_s", "mean"),
            sketch_s=("sketch_mean_s", "mean"),
            speedup=("speedup", "mean"),
            obl_diff=("obl_ratio_diff", "mean"),
        )
        .reset_index()
    )

    compact = compact.rename(columns={
        "edge_bucket": "Number of edges",
        "naive_s": "Naive [s]",
        "sketch_s": "Sketching [s]",
        "speedup": "Speedup",
        "obl_diff": r"$\Delta$ oblivious ratio",
    })

    compact.to_csv(out_dir / "electrical_flow_sketching_summary_compact.csv", index=False)

    # formatted strings for LaTeX
    display_df = compact.copy()
    display_df["Naive [s]"] = display_df["Naive [s]"].map(lambda x: f"{x:.3f}")
    display_df["Sketching [s]"] = display_df["Sketching [s]"].map(lambda x: f"{x:.3f}")
    display_df["Speedup"] = display_df["Speedup"].map(lambda x: f"{x:.2f}x")
    display_df[r"$\Delta$ oblivious ratio"] = display_df[r"$\Delta$ oblivious ratio"].map(
        lambda x: "-" if pd.isna(x) else f"{x:+.3f}"
    )

    # manual LaTeX generation in your preferred style
    lines = []
    lines.append(r"\begin{table}[H]")
    lines.append(r"\centering")
    lines.append(r"\small")
    lines.append(r"\renewcommand{\arraystretch}{1.15}")
    lines.append(r"\begin{tabular}{l|r|r|r|r}")
    lines.append(r"\hline")
    lines.append(r"Number of edges & Naive [s] & Sketching [s] & Speedup & $\Delta$ oblivious ratio \\")
    lines.append(r"\hline")

    for _, row in display_df.iterrows():
        lines.append(
            f"{row['Number of edges']} & "
            f"{row['Naive [s]']} & "
            f"{row['Sketching [s]']} & "
            f"{row['Speedup']} & "
            f"{row[r'$\\Delta$ oblivious ratio']} \\\\"
        )

    lines.append(r"\hline")
    lines.append(r"\end{tabular}")
    lines.append(r"\caption{Summary of the comparison between na\"ive and sketching-based electrical flow algorithms.}")
    lines.append(r"\label{tab:electrical_flow_sketching_summary}")
    lines.append(r"\end{table}")

    with open(out_dir / "electrical_flow_sketching_summary_compact.tex", "w") as f:
        f.write("\n".join(lines) + "\n")

    print("\nCompact sketching vs naive summary:")
    print(display_df.to_string(index=False))


def aggregate_mean_std(df: pd.DataFrame, ycol: str) -> pd.DataFrame:
    g = (
        df.groupby(["solver", "num_edges"], as_index=False)[ycol]
        .agg(["mean", "std", "count"])
        .reset_index()
    )
    g.rename(columns={"count": "n"}, inplace=True)
    g["std"] = g["std"].fillna(0.0)
    return g

def plot_stacked_time_breakdown(
        df: pd.DataFrame,
        solvers: list[str],
        colors: dict[str, str],
        figsize,
        outpath: Path,
):
    """
    Stacked bar plot showing time breakdown by solver (as percentages):
    - solve_time (blue)
    - transformation_time (orange)
    - cycle_removal_time (green)

    Each bar sums to 100%.
    """
    fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)

    # Filter to solvers with data and aggregate
    df_filtered = df[df["solver"].isin(solvers)].copy()

    # Required time components
    time_cols = ["solve_time_micro_seconds", "transformation_time_micro_seconds", "cycle_removal_time_micro_seconds"]

    # Check which columns exist
    available_cols = [col for col in time_cols if col in df_filtered.columns]

    if not available_cols:
        print(f"WARNING: None of the time breakdown columns found. Expected: {time_cols}")
        plt.close(fig)
        return

    # Aggregate by solver
    agg_dict = {col: "mean" for col in available_cols}
    agg_data = df_filtered.groupby("solver", as_index=False)[available_cols].agg(agg_dict)

    # Filter to solvers that have data
    agg_data = agg_data[agg_data["solver"].isin(solvers)]
    if agg_data.empty:
        print("WARNING: No data to plot for time breakdown")
        plt.close(fig)
        return

    # Sort solvers in desired order
    solver_order = [s for s in solvers if s in agg_data["solver"].values]
    agg_data["solver"] = pd.Categorical(agg_data["solver"], categories=solver_order, ordered=True)
    agg_data = agg_data.sort_values("solver")

    # Calculate total time for each solver
    agg_data["total"] = agg_data[available_cols].sum(axis=1)

    # Convert to percentages
    for col in available_cols:
        agg_data[f"{col}_pct"] = (agg_data[col] / agg_data["total"] * 100).fillna(0)

    x = np.arange(len(solver_order))
    width = 0.6

    # Define colors for time components
    component_colors = {
        "solve_time_micro_seconds": "#0072B2",           # blue
        "transformation_time_micro_seconds": "#D55E00",  # orange
        "cycle_removal_time_micro_seconds": "#009E73",   # green
    }

    # Create nicer labels
    component_labels = {
        "solve_time_micro_seconds": "Solve time",
        "transformation_time_micro_seconds": "Transformation time",
        "cycle_removal_time_micro_seconds": "Cycle removal time",
    }

    bottom = np.zeros(len(solver_order))

    for col in available_cols:
        pct_col = f"{col}_pct"
        values = agg_data[pct_col].values
        ax.bar(x, values, width, label=component_labels.get(col, col),
               bottom=bottom, color=component_colors.get(col, "#999999"), alpha=0.8)

        # Add percentage labels on bars
        for i, (val, b) in enumerate(zip(values, bottom)):
            if val > 3:  # Only show label if segment is large enough
                ax.text(x[i], b + val/2, f"{val:.1f}%", ha="center", va="center",
                       fontsize=6, color="white", weight="bold")

        bottom += values

    #ax.set_xlabel("Solver")
    #ax.set_ylabel("Percentage (%)")
    ax.set_ylim(0, 100)
    ax.set_xticks(x)
    #ax.set_xticklabels([pretty_solver_name(s) for s in solver_order], rotation=45, ha="right")
    #ax.legend(loc="upper right", frameon=False, fontsize=7)
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
    fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)

    for s in solvers:
        sub = df[df["solver"] == s].copy()
        if sub.empty:
            continue

        x = pd.to_numeric(sub[xcol], errors="coerce").to_numpy(dtype=float)
        y = pd.to_numeric(sub[ycol], errors="coerce").to_numpy(dtype=float)

        mask = np.isfinite(x) & np.isfinite(y)
        if xlog:
            mask &= (x > 0)
        if ylog:
            mask &= (y > 0)

        x = x[mask]
        y = y[mask]

        if len(x) == 0:
            continue

        if ylog:
            y = np.maximum(y, LOG_EPS)

        # points should dominate
        ax.scatter(
            x, y,
            marker=markers[s],
            s=12,
            alpha=0.7,
            color=colors[s],
            edgecolors="white",
            linewidths=0.2,
            zorder=3,
        )

        # trend line should be subtle
        if add_scaling_line and len(x) >= 3:
            a, b = _fit_powerlaw_line(x, y)
            if a is not None:
                x_min, x_max = np.nanmin(x), np.nanmax(x)
                if xlog:
                    xs = np.logspace(np.log10(x_min), np.log10(x_max), 300)
                else:
                    xs = np.linspace(x_min, x_max, 300)

                ys = a * (xs ** b)
                if ylog:
                    ys = np.maximum(ys, LOG_EPS)

                ax.plot(
                    xs, ys,
                    linestyle="--",
                    linewidth=0.55,
                    color=colors[s],
                    alpha=0.45,
                    zorder=2,
                )

    if xlog:
        ax.set_xscale("log")
    if ylog:
        ax.set_yscale("log")

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.minorticks_on()

    if xlog:
        ax.xaxis.set_major_locator(mticker.LogLocator(base=10))
    else:
        ax.xaxis.set_major_locator(mticker.MaxNLocator(nbins=6))

    if ylog:
        ax.yaxis.set_major_locator(mticker.LogLocator(base=10))
    else:
        ax.yaxis.set_major_locator(mticker.MaxNLocator(nbins=6))

    savefig_all(fig, outpath)
    plt.close(fig)

def parse_arguments():
    """
    Parse command-line arguments for input CSV and output plot directory.
    """
    parser = argparse.ArgumentParser(
        description="Generate publication-grade plots for oblivious routing experiments.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Use default paths
  python3 plot_experiments_pro.py
  
  # Specify custom input CSV and output directory
  python3 plot_experiments_pro.py \\
    --input results/synth/fatclique/combined.csv \\
    --output plots/synthetic/fatclique/
  
  # Specify only input (use default output)
  python3 plot_experiments_pro.py --input my_results.csv
        """
    )

    parser.add_argument(
        "-i", "--input",
        type=str
    )

    parser.add_argument(
        "-o", "--output",
        type=str
    )

    args = parser.parse_args()
    return Path(args.input), Path(args.output)

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

def export_electrical_sketching_overall_table(
        df: pd.DataFrame,
        out_dir: Path,
):
    naive_df = df[df["solver"] == "Electrical Flow (naive)"].copy()
    sketch_df = df[df["solver"] == "Electrical Flow (sketching)"].copy()

    if naive_df.empty or sketch_df.empty:
        return

    naive_mean = naive_df["total_time_micro_seconds"].mean() / 1e6
    sketch_mean = sketch_df["total_time_micro_seconds"].mean() / 1e6

    naive_grouped = naive_df.groupby("num_edges")["total_time_micro_seconds"].mean()
    sketch_grouped = sketch_df.groupby("num_edges")["total_time_micro_seconds"].mean()

    merged = pd.DataFrame({
        "naive": naive_grouped,
        "sketch": sketch_grouped,
    }).dropna()

    merged["speedup"] = merged["naive"] / merged["sketch"]

    overall = pd.DataFrame([{
        "Naive overall mean [s]": f"{naive_mean:.3f}",
        "Sketching overall mean [s]": f"{sketch_mean:.3f}",
        "Overall mean speedup": f"{(naive_mean / sketch_mean):.2f}x",
        "Median grouped speedup": f"{merged['speedup'].median():.2f}x",
        "Best grouped speedup": f"{merged['speedup'].max():.2f}x",
        "Worst grouped speedup": f"{merged['speedup'].min():.2f}x",
    }])

    overall.to_csv(out_dir / "electrical_flow_sketching_overall_summary.csv", index=False)

    with open(out_dir / "electrical_flow_sketching_overall_summary.tex", "w") as f:
        f.write(overall.to_latex(index=False, escape=False))

    print("\nOverall sketching vs naive summary:")
    print(overall.to_string(index=False))

def main():
    set_paper_style()
    RESULT_CSV, OUT_DIR = parse_arguments()
    df = pd.read_csv(RESULT_CSV)
    OUT_DIR = OUT_DIR / "sketching"


    # Clean up graph names: remove dataset prefix (e.g., "Rocketfuel_Topologies/3967.lgf" → "3967")
    if "graph" in df.columns:
        df["graph"] = df["graph"].apply(_graph_short_name)

    if "status" in df.columns:
        # Only filter by status if there are non-NaN values
        if df["status"].notna().any():
            df = df[df["status"] == "OK"].copy()

    all_solvers = sorted(
        [s for s in df["solver"].unique() if "pointer" not in s.lower() or "HST)" in s],
        key=solver_sort_key,
    )

    colors     = solver_palette_unique(all_solvers)
    markers    = solver_markers(all_solvers)
    linestyles = solver_linestyles(all_solvers)

    solvers_electrical = [s for s in all_solvers
                         if "Mendel" not in s and "MendelScaling" not in s
                         and "Raecke" not in s
                          and "Random" not in s
                         and "(Pointer HST)" not in s
                         if "LP" not in s and "Applegate and Cohen" not in s]




    plot_legend_plate(
        solvers_electrical, colors, markers, linestyles,
        outpath=OUT_DIR / "legend_plate",
    )



    agg_runtime = aggregate_mean_std(df[df["solver"].isin(solvers_electrical)].copy(), "total_time_micro_seconds")
    plot_lines(
        agg_runtime, solvers_electrical, colors, markers,
        xcol="num_edges", y_mean_col="mean", y_std_col="std",
        xlabel="Number of edges",
        ylabel="Total running time [microseconds]",
        xlog=True, ylog=True,
        figsize=FIGSIZE_SINGLE,
        outpath=OUT_DIR / "runtime_lines_vs_edges",
        linestyles=linestyles,
    )

    agg_transformation_time = aggregate_mean_std(df[df["solver"].isin(solvers_electrical)].copy(), "transformation_time_micro_seconds")
    plot_lines(
        agg_transformation_time, solvers_electrical, colors, markers,
        xcol="num_edges", y_mean_col="mean", y_std_col="std",
        xlabel="Number of edges",
        ylabel="Tranformation time [microseconds]",
        xlog=False, ylog=True,
        figsize=FIGSIZE_SINGLE,
        outpath=OUT_DIR / "transformation_time_lines_vs_edges",
        linestyles=linestyles,
    )

    agg_mwu_update_time = aggregate_mean_std(df[df["solver"].isin(solvers_electrical)].copy(), "mwu_weight_update_time_micro_seconds")
    plot_lines(
        agg_mwu_update_time, solvers_electrical, colors, markers,
        xcol="num_edges", y_mean_col="mean", y_std_col="std",
        xlabel="Number of edges",
        ylabel="Weight update running time [microseconds]",
        xlog=False, ylog=True,
        figsize=FIGSIZE_SINGLE,
        outpath=OUT_DIR / "weight_update_vs_edges",
        linestyles=linestyles,
    )



    agg_oracle = aggregate_mean_std(df[df["solver"].isin(solvers_electrical)].copy(), "avg_oracle_time_micro_seconds")
    plot_lines(
        agg_oracle, solvers_electrical, colors, markers,
        xcol="num_edges", y_mean_col="mean", y_std_col="std",
        xlabel="Number of edges",
        ylabel="Average oracle running time [microseconds]",
        xlog=False, ylog=True,
        figsize=FIGSIZE_SINGLE,
        outpath=OUT_DIR / "oracle_time_lines_vs_edges",
        linestyles=linestyles,
    )



    df_oblivious = df[df["solver"].isin(solvers_electrical)].dropna(
        subset=["oblivious_ratio"]
    ).copy()
    df_oblivious["oblivious_ratio"] = pd.to_numeric(
        df_oblivious["oblivious_ratio"], errors="coerce"
    )
    df_oblivious = df_oblivious.dropna(subset=["oblivious_ratio"])

    dedup_cols = [c for c in ["graph", "solver","num_nodes", "num_edges", "oblivious_ratio"]
                  if c in df_oblivious.columns]
    df_oblivious_dedup = df_oblivious[dedup_cols].drop_duplicates(
        subset=["graph", "solver"] if "graph" in dedup_cols else None
    )

    solvers_oblivious = [s for s in solvers_electrical
                         if not df_oblivious_dedup[
            df_oblivious_dedup["solver"] == s].empty]

    plot_scatter_cloud(
        df_oblivious_dedup, solvers_oblivious, colors, markers,
        xcol="num_nodes", ycol="oblivious_ratio",
        xlabel="Number of nodes",
        ylabel="Oblivious ratio",
        xlog=False, ylog=False,
        figsize=FIGSIZE_SINGLE,
        outpath=OUT_DIR / "oblivious_ratio_scatter_vs_nodes",
        add_scaling_line=True,
        linestyles=linestyles,
    )

    # Plot electrical flow sketching speedup
    plot_electrical_sketching_speedup(
        df,
        out_dir=OUT_DIR,
        figsize=FIGSIZE_SINGLE,
    )
    plot_electrical_sketching_quality(
        df,
        out_dir=OUT_DIR,
        figsize=FIGSIZE_SINGLE,
    )

    plot_electrical_sketching_quality_gap(
        df,
        out_dir=OUT_DIR,
        figsize=FIGSIZE_SINGLE,
    )

def plot_electrical_sketching_quality_gap(
        df: pd.DataFrame,
        out_dir: Path,
        figsize,
):
    naive_df = df[df["solver"] == "Electrical Flow (naive)"].copy()
    sketch_df = df[df["solver"] == "Electrical Flow (sketching)"].copy()

    for frame in (naive_df, sketch_df):
        frame["num_nodes"] = pd.to_numeric(frame["num_nodes"], errors="coerce")
        frame["oblivious_ratio"] = pd.to_numeric(frame["oblivious_ratio"], errors="coerce")

    naive_grouped = (
        naive_df.dropna(subset=["num_nodes", "oblivious_ratio"])
        .groupby("num_nodes")["oblivious_ratio"]
        .mean()
        .reset_index(name="naive_oblivious_ratio")
    )

    sketch_grouped = (
        sketch_df.dropna(subset=["num_nodes", "oblivious_ratio"])
        .groupby("num_nodes")["oblivious_ratio"]
        .mean()
        .reset_index(name="sketch_oblivious_ratio")
    )

    merged = naive_grouped.merge(sketch_grouped, on="num_nodes", how="inner")
    merged["quality_gap"] = ((merged["sketch_oblivious_ratio"]-merged["naive_oblivious_ratio"])/merged["naive_oblivious_ratio"])*100

    x = merged["num_nodes"].to_numpy(dtype=float)
    y = merged["quality_gap"].to_numpy(dtype=float)

    mask = np.isfinite(x) & np.isfinite(y) & (x > 0)
    x = x[mask]
    y = y[mask]

    fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)

    ax.scatter(
        x, y,
        s=18,
        color="#AA4499",
        edgecolors="white",
        marker=".",
        linewidths=0.3,
        alpha=0.85,
        zorder=3,
    )

    if len(x) >= 2:
        lx = np.log10(x)
        m, b = np.polyfit(lx, y, deg=1)
        xs = np.logspace(np.log10(np.min(x)), np.log10(np.max(x)), 300)
        ys = m * np.log10(xs) + b

        ax.plot(xs, ys, linestyle="--", linewidth=0.9, color="#AA4499", alpha=0.9, zorder=2)

    ax.axhline(0.0, color="gray", linestyle=":", linewidth=0.7, alpha=0.8, zorder=1)

    ax.set_xscale("log")
    #ax.set_xlabel("Number of edges")
    #ax.set_ylabel(r"$\Delta$ oblivious ratio (Sketching - Naive)")
    ax.minorticks_on()
    ax.xaxis.set_major_locator(mticker.LogLocator(base=10))

    savefig_all(fig, out_dir / "electrical_flow_quality_gap_vs_edges")
    plt.close(fig)

def plot_electrical_sketching_speedup(
        df: pd.DataFrame,
        out_dir: Path,
        figsize,
):
    """
    Plot speedup of sketching vs naive for Electrical Flow solvers.
    Speedup = time_naive / time_sketching

    Uses:
    - scatter points (no connecting polyline)
    - best-fit line in log(num_edges)
    """
    naive_df = df[df["solver"] == "Electrical Flow (naive)"].copy()
    sketching_df = df[df["solver"] == "Electrical Flow (sketching)"].copy()

    if naive_df.empty or sketching_df.empty:
        print("WARNING: Missing naive or sketching electrical flow data")
        return

    # Make numeric
    for frame in (naive_df, sketching_df):
        frame["num_edges"] = pd.to_numeric(frame["num_edges"], errors="coerce")
        frame["total_time_micro_seconds"] = pd.to_numeric(
            frame["total_time_micro_seconds"], errors="coerce"
        )

    # Aggregate by num_edges
    naive_grouped = (
        naive_df.groupby("num_edges")["total_time_micro_seconds"]
        .agg(["mean", "count"])
        .reset_index()
    )
    sketching_grouped = (
        sketching_df.groupby("num_edges")["total_time_micro_seconds"]
        .agg(["mean", "count"])
        .reset_index()
    )

    # Merge on num_edges
    speedup_df = naive_grouped.merge(
        sketching_grouped,
        on="num_edges",
        suffixes=("_naive", "_sketching"),
        how="inner",
    )

    if speedup_df.empty:
        print("WARNING: No overlapping num_edges values for naive and sketching")
        return

    speedup_df["speedup"] = speedup_df["mean_naive"] / speedup_df["mean_sketching"]
    speedup_df = speedup_df.sort_values("num_edges")

    x = speedup_df["num_edges"].to_numpy(dtype=float)
    y = speedup_df["speedup"].to_numpy(dtype=float)

    mask = np.isfinite(x) & np.isfinite(y) & (x > 0) & (y > 0)
    x = x[mask]
    y = y[mask]

    if len(x) == 0:
        print("WARNING: No valid speedup points to plot")
        return

    fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)

    # Scatter only
    ax.scatter(
        x, y,
        marker=".",
        s=18,
        color="#0072B2",
        edgecolors="white",
        linewidths=0.3,
        alpha=0.85,
        zorder=3,
    )

    # Best-fit line in log(x)
    if len(x) >= 2:
        lx = np.log10(x)
        m, b = np.polyfit(lx, y, deg=1)

        xs = np.logspace(np.log10(np.min(x)), np.log10(np.max(x)), 300)
        ys = m * np.log10(xs) + b

        ax.plot(
            xs, ys,
            linestyle="--",
            linewidth=0.9,
            color="#D55E00",
            alpha=0.9,
            zorder=2,
        )

        # Optional: print fit info
        print(f"Speedup trend fit: speedup ≈ {m:.3f} * log10(|E|) + {b:.3f}")

    # Reference line: no speedup
    ax.axhline(
        y=1.0,
        color="gray",
        linestyle=":",
        linewidth=0.7,
        alpha=0.8,
        zorder=1,
    )

    ax.set_xscale("log")
    ax.set_yscale("linear")
    #ax.set_xlabel("Number of edges")
    #ax.set_ylabel("Speedup (Naive / Sketching)")
    ax.set_ylim(bottom=0)
    ax.minorticks_on()
    ax.xaxis.set_major_locator(mticker.LogLocator(base=10))

    savefig_all(fig, out_dir / "electrical_flow_sketching_speedup")
    plt.close(fig)

def plot_electrical_sketching_quality(
        df: pd.DataFrame,
        out_dir: Path,
        figsize,
):
    """
    Compare oblivious-ratio quality of Electrical Flow (naive) vs
    Electrical Flow (sketching).

    Uses:
    - grouped scatter points by num_edges
    - best-fit line in log(num_edges) for each variant
    """
    naive_df = df[df["solver"] == "Electrical Flow (naive)"].copy()
    sketch_df = df[df["solver"] == "Electrical Flow (sketching)"].copy()

    if naive_df.empty or sketch_df.empty:
        print("WARNING: Missing naive or sketching electrical flow data")
        return

    for frame in (naive_df, sketch_df):
        frame["num_edges"] = pd.to_numeric(frame["num_edges"], errors="coerce")
        frame["oblivious_ratio"] = pd.to_numeric(frame["oblivious_ratio"], errors="coerce")

    naive_grouped = (
        naive_df.dropna(subset=["num_edges", "oblivious_ratio"])
        .groupby("num_edges")["oblivious_ratio"]
        .agg(["mean", "count"])
        .reset_index()
    )

    sketch_grouped = (
        sketch_df.dropna(subset=["num_edges", "oblivious_ratio"])
        .groupby("num_edges")["oblivious_ratio"]
        .agg(["mean", "count"])
        .reset_index()
    )

    fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)

    series = [
        ("Electrical Flow (naive)", naive_grouped, "#0072B2", "o"),
        ("Electrical Flow (sketching)", sketch_grouped, "#D55E00", "s"),
    ]

    for label, grouped, color, marker in series:
        x = grouped["num_edges"].to_numpy(dtype=float)
        y = grouped["mean"].to_numpy(dtype=float)

        mask = np.isfinite(x) & np.isfinite(y) & (x > 0) & (y > 0)
        x = x[mask]
        y = y[mask]

        if len(x) == 0:
            continue

        ax.scatter(
            x, y,
            marker=marker,
            s=18,
            color=color,
            edgecolors="white",
            linewidths=0.3,
            alpha=0.85,
            zorder=3,
            label=label,
        )

        if len(x) >= 2:
            lx = np.log10(x)
            m, b = np.polyfit(lx, y, deg=1)

            xs = np.logspace(np.log10(np.min(x)), np.log10(np.max(x)), 300)
            ys = m * np.log10(xs) + b

            ax.plot(
                xs, ys,
                linestyle="--",
                linewidth=0.9,
                color=color,
                alpha=0.9,
                zorder=2,
            )

            print(f"{label} quality trend: oblivious_ratio ≈ {m:.3f} * log10(|E|) + {b:.3f}")

    ax.set_xscale("log")
    ax.set_yscale("linear")
    ax.set_xlabel("Number of edges")
    ax.set_ylabel("Oblivious ratio")
    ax.minorticks_on()
    ax.xaxis.set_major_locator(mticker.LogLocator(base=10))

    savefig_all(fig, out_dir / "electrical_flow_quality_vs_edges")
    plt.close(fig)
if __name__ == "__main__":
    main()
