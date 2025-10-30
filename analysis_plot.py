#!/usr/bin/env python3
import argparse
import os
import re
from pathlib import Path
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# ------------------------ Global Settings ------------------------
plt.rcParams.update({
    "font.size": 11,
    "axes.titlesize": 13,
    "axes.labelsize": 11,
    "legend.fontsize": 9,
    "lines.linewidth": 1.4,
    "axes.linewidth": 0.8,
    "legend.frameon": False
})

# Define consistent solver colors
def get_solver_color(name: str) -> str:
    """Return consistent color based on solver name (case-insensitive, robust)."""
    n = name.lower()
    if "electrical" in n:
        return "#1f77b4"  # blue
    elif "tree" in n:
        return "#2ca02c"  # green
    elif "mst" in n:
        return "#7f7f7f"  # dark gray
    elif "cohen" in n:
        return "#8c564b"  # brownish
    elif "random" in n:
        return "#9467bd"  # purple
    elif "lp" in n:
        return "#d62728"  # red
    elif "ckr" in n:
        return "#e377c2"  # pink
    else:
        return "gray"     # fallback

# ------------------------ Helpers ------------------------
def normalize_graph_name(path_str: str) -> str:
    b = os.path.basename(path_str)
    return os.path.splitext(b)[0]

def extract_rep_from_stdout(stdout_path: str):
    m = re.search(r"_rep(\d+)\.out$", os.path.basename(stdout_path))
    return int(m.group(1)) if m else None

def parse_out_metrics(out_path):
    p = Path(out_path)
    try:
        if p.exists():
            p = p.resolve()
    except Exception:
        pass

    run_ms = None
    worst = None
    rand = None
    mwu_iterations = None
    avg_oracle_time = None
    n_nodes = None
    n_edges = None
    demand_congestion = {}  # dict for different demand models
    offline_congestion = None

    # Robust patterns
    re_time = re.compile(r"Running\s*time:\s*([0-9.,]+)\s*\[?milli?seconds\]?", re.IGNORECASE)
    re_worst = re.compile(r"Worst[-\s]*case\s+demand\s+congestion\s*[:=]\s*([0-9.eE+,\-]+)", re.IGNORECASE)
    re_rand_generic = re.compile(r"(random|traffic\s*matrix|tm).*?congestion\s*[:=]\s*([0-9.eE+,\-]+)", re.IGNORECASE)
    re_mwu_iterations = re.compile(r"MWU\s+number\s+of\s+iterations\s*[:=]\s*(\d+)", re.IGNORECASE)
    re_average_oracle_times = re.compile(r"Average\s+oracle\s+time\s*[:=]\s*([0-9.,]+)\s*\[?milli?seconds\]?", re.IGNORECASE)
    re_graph_size = re.compile(r"Graph\s+loaded:\s*(\d+)\s+nodes,\s*(\d+)\s+edges", re.IGNORECASE)
    re_offline_congestion = re.compile(
        r"(?:offline|optimum|optimal)\s+(?:flow|routing|solution|congestion)\s*[:=]\s*([0-9.eE+\-]+)",
        re.IGNORECASE
    )


    # New generic demand line: catch any model name
    re_demand = re.compile(
        r"Worst[-\s]*case\s+congestion\s+for\s+(\w+)\s+model(?:\s+demand)?\s*[:=]\s*([0-9.eE+\-]+)",
        re.IGNORECASE
    )
    #

    try:
        with p.open("r", encoding="utf-8", errors="ignore") as f:
            for raw in f:
                s = raw.strip()

                m = re_graph_size.search(s)
                if m:
                    n_nodes = int(m.group(1))
                    n_edges = int(m.group(2))
                    continue

                m = re_time.search(s)
                if m:
                    run_ms = float(m.group(1).replace(",", "."))
                    continue

                m = re_worst.search(s)
                if m:
                    val = m.group(1).replace(",", ".")
                    try:
                        worst = float(val)
                    except ValueError:
                        worst = None  # or set to a specific error value if needed
                    continue

                m = re_rand_generic.search(s)
                if m:
                    rand = float(m.group(2).replace(",", "."))
                    continue

                m = re_mwu_iterations.search(s)
                if m:
                    mwu_iterations = int(m.group(1))
                    continue

                m = re_average_oracle_times.search(s)
                if m:
                    avg_oracle_time = float(m.group(1).replace(",", "."))
                    continue

                m = re_demand.search(s)
                if m:
                    model = m.group(1).lower()  # e.g. gravity, bimodal, gaussian, uniform
                    val = float(m.group(2).replace(",", "."))
                    demand_congestion[model] = val
                    continue

                m = re_offline_congestion.search(s)
                if m:
                    offline_congestion = float(m.group(1).replace(",", "."))
                    continue




    except FileNotFoundError:
        pass

    return {
        "out_path": str(p),
        "n_nodes": n_nodes,
        "n_edges": n_edges,
        "run_ms": run_ms,
        "worst_congestion": worst,
        "mwu_iterations": mwu_iterations,
        "avg_oracle_time": avg_oracle_time,
        "demand_congestion": demand_congestion,  # nested dict
        "offline_congestion": offline_congestion
    }



def count_nodes_edges_lgf(path):
    """
    Count nodes and edges in .lgf graph
    """
    n, m = None, None
    try:
        with open(path, "r", encoding="utf-8", errors="ignore") as f:
            for line in f:
                s = line.strip()
                if s.startswith("p "):
                    parts = s.split()
                    if len(parts) >= 4:
                        n, m = int(parts[2]), int(parts[3])
                        break
    except Exception:
        pass
    return n, m


# ------------------------ Main ------------------------

def main():
    ap = argparse.ArgumentParser(description="Plots for Oblivious Routing experiments.")
    ap.add_argument("--csv", required=True, help="Path to experiments CSV")
    ap.add_argument("--outdir", default="plots", help="Directory to save plots")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    df = pd.read_csv(args.csv)
    if "graph_name" not in df.columns:
        df["graph_name"] = df["graph"].apply(normalize_graph_name)
    if "rep" not in df.columns and "stdout" in df.columns:
        df["rep"] = df["stdout"].apply(extract_rep_from_stdout)

    # Parse .out logs if needed
    if "stdout" in df.columns:
        df["out_abspath"] = df["stdout"].apply(lambda p: os.path.abspath(p) if isinstance(p, str) and os.path.exists(p) else None)
        metrics = []
        for p in df["out_abspath"].dropna().unique():
            parsed = parse_out_metrics(p)
            if parsed:
                metrics.append(parsed)


        df_out = pd.DataFrame(metrics)
        df = df.merge(df_out, left_on="out_abspath", right_on="out_path", how="left")


        # 🔧 Normalize merged columns
        for col in ["n_nodes", "n_edges", "run_ms", "worst_congestion",
                    "random_max_congestion", "mwu_iterations", "avg_oracle_time"]:
            if f"{col}_x" in df.columns and f"{col}_y" in df.columns:
                df[col] = df[f"{col}_x"].combine_first(df[f"{col}_y"])
                df.drop([f"{col}_x", f"{col}_y"], axis=1, inplace=True, errors="ignore")
            elif f"{col}_x" in df.columns:
                df[col] = df[f"{col}_x"]
                df.drop([f"{col}_x"], axis=1, inplace=True, errors="ignore")
            elif f"{col}_y" in df.columns:
                df[col] = df[f"{col}_y"]
                df.drop([f"{col}_y"], axis=1, inplace=True, errors="ignore")

        print("DEBUG: columns after merge+normalize:", df.columns.tolist())
        print(df[["graph_name", "n_nodes", "n_edges"]].dropna().head())

    # Try to load SciencePlots, fall back to seaborn if missing
    try:
        plt.style.use(['science', 'no-latex'])
    except OSError:
        print("⚠️  'science' style not found — falling back to Seaborn whitegrid.")
        sns.set_theme(style="whitegrid", palette="colorblind")


    # Ensure node/edge counts
    if "n_nodes" not in df.columns or "n_edges" not in df.columns:
        df["graph_lgf"] = df["graph"].apply(lambda p: os.path.abspath(p) if isinstance(p, str) and os.path.exists(p) else None)
        df[["n_nodes", "n_edges"]] = df["graph_lgf"].apply(lambda p: pd.Series(count_nodes_edges_lgf(p) if p else (None, None)))

    # Save merged for debugging
    df.to_csv(os.path.join(args.outdir, "merged_debug.csv"), index=False)


    # --------- Plot 1: Average running time by solver ---------
    if "run_ms" in df.columns:
        agg = df.groupby("solver")["run_ms"].agg(["mean", "std"])
        fig, ax = plt.subplots()
        agg["mean"].div(1000).plot.bar(yerr=agg["std"].div(1000), ax=ax, capsize=4)
        ax.set_title("Average Running Time by Solver")
        ax.set_ylabel("Time (s)")
        fig.tight_layout()
        plt.savefig(os.path.join(args.outdir, "avg_runtime_by_solver.png"), dpi=180)
        plt.close()

    # --------- Plot 2: Runtime scaling by solver (log–log, clean) ---------
    # robust style fallback (no hard dependency on SciencePlots)
    try:
        plt.style.use(['science', 'no-latex'])
    except OSError:
        sns.set_theme(style="whitegrid", palette="colorblind")

    # extend colors to match your solver labels exactly
    colors = {
        "electrical": "#1f77b4",
        "electrical_optimized": "#1f77b4",
        "tree": "#2ca02c",
        "mst": "#7f7f7f",
        "cohen": "#8c8c8c",
        "ckr": "#e377c2",
        "lp": "#d62728",
    }

    if "n_edges" in df.columns and "run_ms" in df.columns:
        runs = df.dropna(subset=["n_edges", "run_ms"]).copy()
        # keep per-instance (graph_name, solver) means to de-noise replicates,
        # but don't average across different graphs
        runs = runs.groupby(["graph_name", "solver"], as_index=False).agg({
            "run_ms": "mean",
            "n_edges": "first"
        })

        fig, ax = plt.subplots(figsize=(6.8, 3.8), dpi=180)

        # plot per-solver
        for solver, sub in runs.groupby("solver"):
            sub = sub.sort_values("n_edges")
            name = solver.lower()
            color = get_solver_color(name)

            # scatter points
            s = 18 if name != "electrical_optimized" and name != "electrical" else 28
            ec = "none" if name != "electrical_optimized" and name != "electrical" else "black"
            lw = 0 if name != "electrical_optimized" and name != "electrical" else 0.4
            z = 2 if name != "electrical_optimized" and name != "electrical" else 4
            alpha = 0.8 if name != "electrical_optimized" and name != "electrical" else 1.0

            ax.scatter(
                sub["n_edges"], sub["run_ms"]/1000.0,
                s=s, color=color, alpha=alpha, zorder=z,
                label=solver if name not in {"electrical", "electrical_optimized"} else f"{solver} (ours)",
                edgecolor=ec, linewidth=lw, marker="o"
            )

            # power-law fit (straight line in log–log)
            if len(sub) >= 3:
                x = sub["n_edges"].to_numpy()
                y = (sub["run_ms"]/1000.0).to_numpy()
                # filter out non-positive (just in case)
                mask = (x > 0) & (y > 0)
                x, y = x[mask], y[mask]
                if len(x) >= 3:
                    a, b = np.polyfit(np.log10(x), np.log10(y), 1)   # log10 y = a log10 x + b
                    xx = np.linspace(x.min(), x.max(), 100)
                    yy = 10**b * (xx**a)
                    ax.plot(xx, yy,
                            color=color,
                            alpha=0.7 if name in {"electrical", "electrical_optimized"} else 0.45,
                            linewidth=1.6 if name in {"electrical", "electrical_optimized"} else 1.0,
                            zorder=1)

                    # show slope for our method
                    if name in {"electrical", "electrical_optimized"}:
                        ax.text(0.04, 0.92, f"Slope ≈ {a:.2f}",
                                transform=ax.transAxes, fontsize=9, color=color, ha="left", va="top")

        # log–log axes
        ax.set_xscale("log")
        ax.set_yscale("log")

        # axis labels / title
        ax.set_title("Runtime Scaling by Solver")
        ax.set_xlabel("Number of edges |E| (log scale)")
        ax.set_ylabel("Runtime (seconds, log scale)")

        # lighter grid: major faint, minor very faint
        ax.grid(True, which="major", linestyle="--", alpha=0.25)
        ax.grid(True, which="minor", linestyle=":", alpha=0.12)

        # legend outside to avoid covering points
        ax.legend(frameon=False, loc="center left", bbox_to_anchor=(1.02, 0.5))
        fig.tight_layout()
        plt.savefig(os.path.join(args.outdir, "runtime_per_instance_thesis.svg"),
                    bbox_inches="tight")
        plt.close()



    # --------- Plot 3: Average oracle running time (Electrical & Tree only) ---------
    if "avg_oracle_time" in df.columns and "n_edges" in df.columns:
        runs = df.dropna(subset=["avg_oracle_time", "n_edges"]).copy()
        # exclude LP-based solvers
        runs = runs.loc[~runs["solver"].str.lower().str.contains("lp")]
        runs = runs.groupby(["graph_name", "solver"], as_index=False).agg({
            "avg_oracle_time": "mean",
            "n_edges": "first"   # preserve edge count for sorting
        })
        runs = runs.sort_values("n_edges")

        fig, ax = plt.subplots()

        for solver, sub in runs.groupby("solver"):
            if solver == "cohen":
                continue  # skip cohen for clarity
            mean = sub["avg_oracle_time"].mean()
            ax.scatter(
                sub["n_edges"],
                sub["avg_oracle_time"],
                s=5,            # 👈 small dot size
                label=f"{solver} avg = {mean:.3f} ms",
                marker="o"       # keep circular dots
            )

        ax.set_title("Average oracle running time per Instance (sorted by |E|)")
        ax.set_ylabel("Running time (ms)")
        # ax.set_xticklabels([])   # hide the text labels
        ax.legend()
        fig.autofmt_xdate(rotation=90)
        fig.tight_layout()
        plt.savefig(os.path.join(args.outdir, "oracle_avg_instance.png"), dpi=180)
        plt.close()


    # --------- Plot 4: MWU iterations per instance (sorted by edges) ---------
    if "mwu_iterations" in df.columns and "n_edges" in df.columns:
        runs = df.dropna(subset=["mwu_iterations", "n_edges"]).copy()
        # exclude LP-based solvers
        runs = runs.loc[~runs["solver"].str.lower().str.contains("lp")]
        runs = runs.groupby(["graph_name", "solver"], as_index=False).agg({
            "mwu_iterations": "mean",
            "n_edges": "first"   # preserve edge count for sorting
        })
        runs = runs.sort_values("n_edges")

        fig, ax = plt.subplots()

        for solver, sub in runs.groupby("solver"):
            if solver == "cohen":
                continue  # skip cohen for clarity
            mean = sub["mwu_iterations"].mean()
            ax.scatter(
                sub["n_edges"],
                sub["mwu_iterations"],
                s=5,            # 👈 small dot size
                label=f"{solver} avg = {mean:.3f}",
                marker="o"       # keep circular dots
            )

        ax.set_title("MWU Iterations per Instance (sorted by |E|)")
        ax.set_ylabel("Iterations")
        # ax.set_xticklabels([])   # hide the text labels
        ax.legend()
        fig.autofmt_xdate(rotation=90)
        fig.tight_layout()
        plt.savefig(os.path.join(args.outdir, "mwu_iterations_per_instance.png"), dpi=180)
        plt.close()





    # --------- Relative congestion (expanded) ---------
    if "demand_congestion" in df.columns:
        demand_df = df["demand_congestion"].apply(lambda x: pd.Series(x) if isinstance(x, dict) else pd.Series({}))
        df_expanded = pd.concat([df, demand_df], axis=1)

        keys = ["graph_name", "rep"] if "rep" in df_expanded.columns else ["graph_name"]



        df_expanded["offline_congestion"] = pd.to_numeric(df_expanded["offline_congestion"], errors="coerce")

        # --- compute relative % columns ---
        # --------- 📊 Relative % Error per Instance (sorted by |E|) ---------
        for demand in ["bimodal", "gaussian", "gravity", "uniform"]:
            if demand not in df_expanded.columns:
                continue
            df_expanded[demand] = pd.to_numeric(df_expanded[demand], errors="coerce")

            err_col = f"{demand}_rel_err_pct"

            num = df_expanded[demand].to_numpy(dtype=float, copy=False)
            den = df_expanded["offline_congestion"].to_numpy(dtype=float, copy=False)

            valid = ~np.isnan(num) & ~np.isnan(den) & (den != 0.0)
            err = np.full(len(df_expanded), np.nan, dtype=float)
            np.divide(np.abs(num - den), np.abs(den), out=err, where=valid)
            df_expanded[err_col] = err

            runs = df_expanded.dropna(subset=[err_col, "solver", "graph_name", "n_edges"]).copy()
            if runs.empty:
                continue

            # Group and sort by number of edges
            runs = (
                runs.groupby(["graph_name", "solver"], as_index=False)
                .agg({err_col: "mean", "n_edges": "first"})
                .sort_values("n_edges")
            )

            # Create scatter plot
            fig, ax = plt.subplots(figsize=(12, 4))
            # Group by solver and by instance

            for solver, sub in runs.groupby("solver"):
                avg_err = sub[err_col].mean()
                ax.scatter(
                    sub["graph_name"],
                    sub[err_col],
                    s=15,
                    label=f"{solver} avg = {avg_err:.3f}%",
                )



            ax.axhline(0, color="black", linestyle="--", linewidth=1)
            ax.set_title(f"{demand.capitalize()} – Relative error from offline optimum (|offline optimal routing - solver routing|/offline optimal routing)")
            ax.set_ylabel("Relative error in %")
            ax.set_xlabel("Instances sorted by |E|")
            ax.set_xticklabels([])   # hide the text labels
            ax.legend()
            fig.tight_layout()
            plt.savefig(os.path.join(args.outdir, f"rel_error_from_offline_{demand}_per_instance.png"), dpi=180)
            plt.close()

        # --------- Plot 5: Quality vs Runtime Tradeoff (Efficiency Frontier) ---------
    # --------- Plot 5: Quality vs Runtime (all instances, not averaged) ---------
    if "run_ms" in df.columns and "demand_congestion" in df.columns:
        # expand demand models
        demand_df = df["demand_congestion"].apply(
            lambda x: pd.Series(x) if isinstance(x, dict) else pd.Series({})
        )
        df_expanded = pd.concat([df, demand_df], axis=1)

        demand_metric = "gravity"  # pick your demand model here
        if demand_metric not in df_expanded.columns:
            print(f"⚠️ Demand model '{demand_metric}' not found, skipping quality-runtime plot.")
        else:
            df_expanded[demand_metric] = pd.to_numeric(df_expanded[demand_metric], errors="coerce")
            df_expanded["offline_congestion"] = pd.to_numeric(df_expanded["offline_congestion"], errors="coerce")

            valid = (~df_expanded["offline_congestion"].isna()) & (df_expanded["offline_congestion"] != 0)
            df_expanded["rel_err_pct"] = np.nan
            df_expanded.loc[valid, "rel_err_pct"] = (
                    np.abs(df_expanded[demand_metric] - df_expanded["offline_congestion"])
                    / df_expanded["offline_congestion"] * 100
            )

            # prepare data
            runs = df_expanded.dropna(subset=["run_ms", "rel_err_pct", "solver"]).copy()
            runs["runtime_s"] = runs["run_ms"] / 1000.0

            # compute averages for highlighting
            agg = (
                runs.groupby("solver")
                .agg(avg_runtime_s=("runtime_s", "mean"),
                     avg_error_pct=("rel_err_pct", "mean"))
                .reset_index()
            )


            sns.set_theme(style="whitegrid", palette="colorblind")

            fig, ax = plt.subplots(figsize=(6.4, 4.2), dpi=180)

            # --- per-instance scatter cloud ---
            for solver, sub in runs.groupby("solver"):
                name = str(solver)
                color = get_solver_color(name)
                ax.scatter(sub["runtime_s"], sub["rel_err_pct"],
                           s=18, alpha=0.4, color=color, label=None, zorder=1)

            # --- average markers (highlighted) ---
            for _, row in agg.iterrows():
                name = row["solver"].lower()
                color = get_solver_color(name)
                edge = "black" if "electrical" in name else "none"
                lw = 0.5 if "electrical" in name else 0
                ax.scatter(row["avg_runtime_s"], row["avg_error_pct"],
                           s=110, color=color, edgecolor=edge, linewidth=lw,
                           zorder=3)
                ax.text(row["avg_runtime_s"] * 1.08,
                        row["avg_error_pct"] * 1.08,
                        row["solver"],
                        fontsize=9, color=color)

            # --- annotate our method ---
            if "electrical" in agg["solver"].str.lower().values:
                elec = agg[agg["solver"].str.lower().str.contains("electrical")].iloc[0]
                ax.annotate("Electrical Flow (ours)",
                            (elec["avg_runtime_s"], elec["avg_error_pct"]),
                            xytext=(25, -15), textcoords="offset points",
                            arrowprops=dict(arrowstyle="->", lw=1, color="#1f77b4"),
                            fontsize=9, color="#1f77b4", weight="bold")

            # --- axes, scales, cosmetics ---
            ax.set_xscale("log")
            ax.set_xlabel("Runtime (seconds, log scale)")
            ax.set_ylabel("Relative error (%)")
            ax.set_title(f"Solution Quality vs Runtime ({demand_metric.capitalize()} demands)")

            ax.grid(True, which="major", linestyle="--", alpha=0.3)
            ax.grid(True, which="minor", linestyle=":", alpha=0.12)
            ax.text(0.96, 0.04, "← Faster / Better →", transform=ax.transAxes,
                    fontsize=9, ha="right", va="bottom", color="gray")

            # optional regression per solver (thin lines)
            for solver, sub in runs.groupby("solver"):
                if len(sub) >= 5:
                    color = colors.get(solver.lower(), "gray")
                    x = sub["runtime_s"].to_numpy()
                    y = sub["rel_err_pct"].to_numpy()
                    mask = (x > 0) & (y > 0)
                    x, y = x[mask], y[mask]
                    if len(x) >= 3:
                        a, b = np.polyfit(np.log10(x), np.log10(y), 1)
                        xx = np.linspace(x.min(), x.max(), 100)
                        yy = 10**b * (xx**a)
                        ax.plot(xx, yy, color=color, alpha=0.4, linewidth=1.0, zorder=0)

            fig.tight_layout()
            plt.savefig(os.path.join(args.outdir, f"quality_vs_runtime_all_{demand_metric}.svg"),
                        bbox_inches="tight")
            plt.close()






if __name__ == "__main__":
    main()
