#!/usr/bin/env python3
import os, re, argparse
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

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

    # Robust patterns
    re_time = re.compile(r"Running\s*time:\s*([0-9.,]+)\s*\[?milli?seconds\]?", re.IGNORECASE)
    re_worst = re.compile(r"Worst[-\s]*case\s+demand\s+congestion\s*[:=]\s*([0-9.eE+,\-]+)", re.IGNORECASE)
    re_rand_generic = re.compile(r"(random|traffic\s*matrix|tm).*?congestion\s*[:=]\s*([0-9.eE+,\-]+)", re.IGNORECASE)
    re_mwu_iterations = re.compile(r"MWU\s+number\s+of\s+iterations\s*[:=]\s*(\d+)", re.IGNORECASE)
    re_average_oracle_times = re.compile(r"Average\s+oracle\s+time\s*[:=]\s*([0-9.,]+)\s*\[?milli?seconds\]?", re.IGNORECASE)
    re_graph_size = re.compile(r"Graph\s+loaded:\s*(\d+)\s+nodes,\s*(\d+)\s+edges", re.IGNORECASE)

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
                    worst = float(m.group(1).replace(",", "."))
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
        "demand_congestion": demand_congestion  # nested dict
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

    # --------- Plot 2: Running time per instance sorted by nodes ---------
    if "n_edges" in df.columns and "run_ms" in df.columns:
        runs = df.dropna(subset=["n_edges", "run_ms"]).copy()
        runs = runs.groupby(["graph_name", "solver"], as_index=False).agg({
            "run_ms": "mean",
            "n_edges": "first"   # keep node count for sorting
        })

        runs = runs.sort_values("n_edges")
        fig, ax = plt.subplots()
        for solver, sub in runs.groupby("solver"):
            ax.scatter(
                sub["n_edges"],
                sub["run_ms"],  # convert ms to s
                s=5,            # 👈 small dot size
                label=solver,
                marker="o"       # keep circular dots
            )

        ax.set_title("Running Time per Instance (sorted by |E|)")
        ax.set_ylabel("Time (ms)")
        # ax.set_xticklabels([])   # hide the text labels
        ax.legend()
        fig.autofmt_xdate(rotation=90)
        fig.tight_layout()
        plt.savefig(os.path.join(args.outdir, "runtime_per_instance.png"), dpi=180)
        plt.close()

    # --------- Plot 3: Average oracle running time (Electrical & Tree only) ---------
    if "avg_oracle_time" in df.columns:
        oracles = df.dropna(subset=["avg_oracle_time"])
        agg = oracles.groupby("solver")["avg_oracle_time"].agg(["mean", "std"])
        fig, ax = plt.subplots()
        agg["mean"].div(1000).plot.bar(yerr=agg["std"].div(1000), ax=ax, capsize=4)
        ax.set_title("Average Oracle Time ")
        ax.set_ylabel("Time (s)")
        fig.tight_layout()
        plt.savefig(os.path.join(args.outdir, "avg_oracle_time.png"), dpi=180)
        plt.close()

    # --------- Plot 4: MWU iterations per instance (sorted by edges) ---------
    if "mwu_iterations" in df.columns and "n_edges" in df.columns:
        runs = df.dropna(subset=["mwu_iterations", "n_edges"]).copy()
        runs = runs.groupby(["graph_name", "solver"], as_index=False).agg({
            "mwu_iterations": "mean",
            "n_edges": "first"   # preserve edge count for sorting
        })
        runs = runs.sort_values("n_edges")

        fig, ax = plt.subplots()

        for solver, sub in runs.groupby("solver"):
            ax.scatter(
                sub["n_edges"],
                sub["mwu_iterations"],
                s=5,            # 👈 small dot size
                label=solver,
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

    # --------- Plot 5: Worst-case congestion per demand model ---------
    if "demand_congestion" in df.columns:
        # Expand dict column into multiple demand columns
        demand_df = df["demand_congestion"].dropna().apply(pd.Series)
        df = pd.concat([df, demand_df], axis=1)

        for demand in ["bimodal", "gaussian", "gravity", "uniform"]:
            if demand in df.columns:
                runs = df.dropna(subset=[demand]).copy()
                runs = runs.groupby(["graph_name", "solver"])[demand].mean().reset_index()
                fig, ax = plt.subplots()
                for solver, sub in runs.groupby("solver"):
                    ax.plot(sub["graph_name"], sub[demand], label=solver)
                ax.set_title(f"Worst-case Congestion ({demand})")
                ax.set_ylabel("Congestion")
                ax.set_xticklabels([])   # hide the text labels
                ax.legend()
                fig.autofmt_xdate(rotation=90)
                fig.tight_layout()
                plt.savefig(os.path.join(args.outdir, f"worst_congestion_{demand}.png"), dpi=180)
                plt.close()

        # --------- Extra Plot: Average worst-case congestion per demand model ---------
    if "demand_congestion" in df.columns:
        demand_df = df["demand_congestion"].dropna().apply(
            lambda x: pd.Series(x) if isinstance(x, dict) else pd.Series({})
        )
        df_expanded = pd.concat([df, demand_df], axis=1)

        averages = {}
        for demand in ["bimodal", "gaussian", "gravity", "uniform"]:
            if demand in df_expanded.columns:
                # Extract only numeric values safely
                vals = df_expanded[demand].apply(
                    lambda x: float(x) if isinstance(x, (int, float)) else None
                ).dropna()
                if not vals.empty:
                    averages[demand] = vals.mean()

        if averages:
            fig, ax = plt.subplots()
            pd.Series(averages, dtype=float).plot.bar(ax=ax, color="skyblue")
            ax.set_title("Average Worst-case Congestion by Demand Model")
            ax.set_ylabel("Congestion")
            ax.set_xlabel("Demand Model")
            fig.tight_layout()
            plt.savefig(os.path.join(args.outdir, "avg_worst_congestion_per_demand.png"), dpi=180)
            plt.close()


if __name__ == "__main__":
    main()
