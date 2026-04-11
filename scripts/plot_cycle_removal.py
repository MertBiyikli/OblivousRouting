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



    if xlog:
        ax.set_xscale("log")
    if ylog:
        ax.set_yscale("log")

    # _apply_labels(ax, xlabel, ylabel)
    ax.minorticks_on()

    # Show more x-axis ticks for better readability without clutter
    if xlog:
        # For log scale, use LogLocator to respect logarithmic spacing
        ax.xaxis.set_major_locator(mticker.LogLocator(base=10, numticks=6))
    else:
        # For linear scale, use MaxNLocator with moderate number of ticks
        ax.xaxis.set_major_locator(mticker.MaxNLocator(nbins=6, integer=False))



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
                    xs = np.logspace(np.log10(x_min), np.log10(x_max), 500)
                else:
                    xs = np.linspace(x_min, x_max, 500)
                ys = a * (xs ** b)
                if ylog:
                    ys = np.maximum(ys, LOG_EPS)
                ax.plot(xs, ys, linestyle="-", linewidth=0.7, color=colors[s], alpha=0.85, zorder=2)

    if xlog:
        ax.set_xscale("log")
    if ylog:
        ax.set_yscale("log")

    ax.minorticks_on()


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

def plot_cycle_removal_comparison(
        df: pd.DataFrame,
        solvers: list[str],
        colors: dict[str, str],
        markers: dict[str, str],
        out_dir: Path,
        figsize,
        linestyles: dict[str, str] | None = None,
):
    """
    For each tree solver, create a line plot comparing runtime between
    different cycle removal strategies.
    Only plots for Raecke FRT, Raecke CKR, and Random MST solvers and MendelScaling variants.
    """
    # Check if cycle_removal_type column exists
    if "cycle_removal_type" not in df.columns:
        print("WARNING: cycle_removal_type column not found in data")
        return

    # Filter to only Raecke FRT, CKR, and MST solvers with Flat HST
    tree_solver_keywords = ["Raecke FRT", "Raecke CKR", "Random MST","Raecke FRT + MendelScaling", "Raecke CKR + MendelScaling" ]
    filtered_solvers = [s for s in solvers if any(keyword in s for keyword in tree_solver_keywords) and "(Flat HST)" in s]

    # Get unique solvers (base names)
    unique_solvers = set(filtered_solvers)
    print(f"Plotting cycle removal comparison for: {unique_solvers}")

    for solver in sorted(unique_solvers):
        # Get all data for this solver with all strategies
        solver_df = df[df["solver"] == solver].copy()
        if solver_df.empty:
            continue

        # Get unique strategies for this solver
        strategies = sorted(solver_df["cycle_removal_type"].dropna().unique())
        if len(strategies) < 2:
            continue  # Skip if only one strategy

        fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)

        for strategy in strategies:
            strategy_df = solver_df[solver_df["cycle_removal_type"] == strategy].copy()
            if strategy_df.empty:
                continue

            # Aggregate by num_edges
            agg = strategy_df.groupby("num_edges")["cycle_removal_time_micro_seconds"].agg(["mean", "std"]).reset_index()
            agg.columns = ["num_edges", "mean", "std"]
            agg = agg.sort_values("num_edges")

            x = agg["num_edges"].to_numpy()
            y = agg["mean"].to_numpy()

            # Apply log scale for y values
            y = np.maximum(y, LOG_EPS)

            # Create a unique color and marker for each strategy
            strategy_idx = list(strategies).index(strategy)
            color = _PALETTE[strategy_idx % len(_PALETTE)]
            marker = solver_markers([f"dummy_{i}" for i in range(len(strategies))])[f"dummy_{strategy_idx}"]
            ls = _LINESTYLES[strategy_idx % len(_LINESTYLES)]

            ax.plot(
                x, y,
                linestyle=ls,
                marker=marker,
                markersize=2.5,
                markeredgewidth=0.4,
                markeredgecolor="white",
                label=strategy,
                color=color,
                linewidth=0.9,
            )

        ax.set_xscale("linear")
        ax.set_yscale("log")
        #ax.set_xlabel("Number of edges")
        #ax.set_ylabel("Cycle removal time [microseconds]")
        ax.minorticks_on()
        # ax.legend(loc="best", frameon=False, fontsize=8)

        # Use solver name for filename
        solver_name = solver.replace(" ", "_").replace("(", "").replace(")", "").lower()
        savefig_all(fig, out_dir / f"cycle_removal_speedup_{solver_name}")
        plt.close(fig)

def plot_cycle_removal_strategy_speedup(
        df: pd.DataFrame,
        out_dir: Path,
        figsize,
):
    """
    Plot speedup of Tarjan_SCC vs naive cycle removal strategies across all solvers.
    Speedup = time_naive / time_tarjan
    Combines all tree solvers (FRT, CKR, MST and their variants).
    """
    # Check if cycle_removal_type column exists
    if "cycle_removal_type" not in df.columns:
        print("WARNING: cycle_removal_type column not found in data")
        return

    # Filter to tree solvers only
    tree_solver_keywords = ["Raecke FRT", "Raecke CKR", "Random MST"]
    df_tree = df[df["solver"].str.contains("|".join(tree_solver_keywords), na=False)].copy()

    if df_tree.empty:
        print("WARNING: No tree solvers found")
        return

    # Get cycle removal strategies
    strategies = sorted(df_tree["cycle_removal_type"].dropna().unique())

    # Identify naive and tarjan
    naive_df = df_tree[df_tree["cycle_removal_type"] == "naive"].copy()
    tarjan_df = df_tree[df_tree["cycle_removal_type"] == "Tarjan_SCC"].copy()

    if naive_df.empty or tarjan_df.empty:
        print(f"WARNING: Missing naive or Tarjan_SCC data. Available strategies: {strategies}")
        return

    print(f"Computing speedup: naive vs Tarjan_SCC across {len(naive_df)} + {len(tarjan_df)} instances")

    # Group by num_edges and compute average times
    naive_grouped = naive_df.groupby("num_edges")["total_time_micro_seconds"].agg(["mean", "count"]).reset_index()
    tarjan_grouped = tarjan_df.groupby("num_edges")["total_time_micro_seconds"].agg(["mean", "count"]).reset_index()

    # Merge on num_edges to align instances
    speedup_df = naive_grouped.merge(
        tarjan_grouped, on="num_edges", suffixes=("_naive", "_tarjan")
    )

    # Compute speedup: speedup = time_naive / time_tarjan
    speedup_df["speedup"] = speedup_df["mean_naive"] / speedup_df["mean_tarjan"]

    # Sort by num_edges
    speedup_df = speedup_df.sort_values("num_edges")

    x = speedup_df["num_edges"].to_numpy()
    y = speedup_df["speedup"].to_numpy()

    fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)

    ax.plot(
        x, y,
        linestyle="-",
        marker="o",
        markersize=2.5,
        markeredgewidth=0.4,
        markeredgecolor="white",
        label="Tarjan_SCC vs Naive Speedup",
        color="#009E73",
        linewidth=0.9,
    )

    # Add horizontal line at y=1 (no speedup)
    ax.axhline(y=1.0, color="gray", linestyle=":", linewidth=0.7, alpha=0.7, label="No speedup")

    ax.set_xscale("linear")
    ax.set_yscale("log")
    #ax.set_xlabel("Number of edges")
    #ax.set_ylabel("Speedup (Naive / Tarjan_SCC)")
    ax.set_ylim(bottom=0)
    ax.minorticks_on()
    #ax.legend(loc="best", frameon=False, fontsize=8)
    ax.grid(True, alpha=0.3)

    savefig_all(fig, out_dir / "cycle_removal_strategy_speedup")
    plt.close(fig)

def plot_electrical_sketching_speedup(
        df: pd.DataFrame,
        out_dir: Path,
        figsize,
):
    """
    Plot speedup of sketching vs naive for Electrical Flow solvers.
    Speedup = time_naive / time_sketching
    """
    # Filter to electrical flow solvers
    electrical_solvers = [s for s in df["solver"].unique()
                          if "Electrical Flow" in s]

    if len(electrical_solvers) < 2:
        print("WARNING: Need at least 2 electrical flow variants for speedup comparison")
        return

    print(f"Computing speedup for: {electrical_solvers}")

    # Identify naive and sketching versions
    naive_df = df[df["solver"] == "Electrical Flow (naive)"].copy()
    sketching_df = df[df["solver"] == "Electrical Flow (sketching)"].copy()

    if naive_df.empty or sketching_df.empty:
        print("WARNING: Missing naive or sketching electrical flow data")
        return

    # Group by num_edges and compute average times
    naive_grouped = naive_df.groupby("num_edges")["total_time_micro_seconds"].agg(["mean", "count"]).reset_index()
    sketching_grouped = sketching_df.groupby("num_edges")["total_time_micro_seconds"].agg(["mean", "count"]).reset_index()

    # Merge on num_edges to align instances
    speedup_df = naive_grouped.merge(
        sketching_grouped, on="num_edges", suffixes=("_naive", "_sketching")
    )

    # Compute speedup: speedup = time_naive / time_sketching
    speedup_df["speedup"] = speedup_df["mean_naive"] / speedup_df["mean_tarjan"]

    # Sort by num_edges
    speedup_df = speedup_df.sort_values("num_edges")

    x = speedup_df["num_edges"].to_numpy()
    y = speedup_df["speedup"].to_numpy()

    fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)

    ax.plot(
        x, y,
        linestyle="-",
        marker="o",
        markersize=2.5,
        markeredgewidth=0.4,
        markeredgecolor="white",
        label="Sketching vs Naive Speedup",
        color="#0072B2",
        linewidth=0.9,
    )

    # Add horizontal line at y=1 (no speedup)
    ax.axhline(y=1.0, color="gray", linestyle=":", linewidth=0.7, alpha=0.7, label="No speedup")

    ax.set_xscale("linear")
    ax.set_yscale("log")
    #ax.set_xlabel("Number of edges")
    #ax.set_ylabel("Speedup (Naive / Sketching)")
    ax.set_ylim(bottom=0)
    ax.minorticks_on()
    #ax.legend(loc="best", frameon=False, fontsize=8)
    ax.grid(True, alpha=0.3)

    savefig_all(fig, out_dir / "cycle_removal_speedup")
    plt.close(fig)

def main():
    set_paper_style()
    RESULT_CSV, OUT_DIR = parse_arguments()
    df = pd.read_csv(RESULT_CSV)
    OUT_DIR = OUT_DIR / "cycle"

    if "graph" in df.columns:
        df["graph"] = df["graph"].apply(_graph_short_name)

    if "status" in df.columns:
        # Only filter by status if there are non-NaN values
        if df["status"].notna().any():
            df = df[df["status"] == "OK"].copy()

    all_solvers = sorted(
        [s for s in df["solver"].unique() if "(Pointer HST)" not in s.lower()],
        key=solver_sort_key,
    )


    colors     = solver_palette_unique(all_solvers)
    markers    = solver_markers(all_solvers)
    linestyles = solver_linestyles(all_solvers)


    print(f"All solvers: {all_solvers}")

    # Plot cycle removal comparison for all solvers
    plot_cycle_removal_comparison(
        df, all_solvers, colors, markers,
        out_dir=OUT_DIR,
        figsize=FIGSIZE_SINGLE,
        linestyles=linestyles,
    )

    # Plot cycle removal strategy speedup (naive vs Tarjan) across all tree solvers
    plot_cycle_removal_strategy_speedup(
        df,
        out_dir=OUT_DIR,
        figsize=FIGSIZE_SINGLE,
    )



if __name__ == "__main__":
    main()
