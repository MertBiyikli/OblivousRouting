#!/usr/bin/env python3
import os, re, argparse
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

# ------------------------ Helpers ------------------------

def find_solver_col(columns, key="cohen"):
    """Return the column in `columns` whose name contains `key` (case-insensitive), else None."""
    for c in columns:
        if key in str(c).lower():
            return c
    return None


def normalize_graph_name(path_str: str) -> str:
    b = os.path.basename(path_str)
    return os.path.splitext(b)[0]

def extract_rep_from_stdout(stdout_path: str):
    m = re.search(r"_rep(\d+)\.out$", os.path.basename(stdout_path))
    return int(m.group(1)) if m else None

def parse_out_metrics(out_path):
    """
    Accepts str or Path. Returns dict with out_path (absolute),
    run_ms, worst_congestion, random_max_congestion.
    Recognizes:
      - 'Running time: 6067 [milliseconds]'
      - 'Worst case demand congestion: 228'
      - 'Worst case congestion for gravity model demand: 37.2502'
      - generic '(random|traffic matrix|tm) ... congestion: X'
    """
    p = Path(out_path)
    try:
        if p.exists():
            p = p.resolve()
    except Exception:
        pass

    run_ms = None
    worst = None
    rand = None

    # Robust patterns
    re_time = re.compile(r"Running\s*time:\s*([0-9.,]+)\s*\[?milli?seconds\]?", re.IGNORECASE)
    re_worst = re.compile(r"Worst[-\s]*case\s+demand\s+congestion\s*[:=]\s*([0-9.eE+,\-]+)", re.IGNORECASE)
    re_gravity_exact = re.compile(
        r"Worst[-\s]*case\s+congestion\s+for\s+gravity\s+model(?:\s+demand)?\s*[:=]\s*([0-9.eE+,\-]+)",
        re.IGNORECASE
    )
    re_rand_generic = re.compile(
        r"(random|traffic\s*matrix|tm).*?congestion\s*[:=]\s*([0-9.eE+,\-]+)", re.IGNORECASE
    )

    try:
        with p.open("r", encoding="utf-8", errors="ignore") as f:
            for raw in f:
                s = raw.strip()

                m = re_time.search(s)
                if m:
                    run_ms = float(m.group(1).replace(",", "."))
                    continue

                m = re_worst.search(s)
                if m:
                    worst = float(m.group(1).replace(",", "."))
                    continue

                m = re_gravity_exact.search(s)
                if m:
                    rand = float(m.group(1).replace(",", "."))
                    continue

                m = re_rand_generic.search(s)
                if m:
                    rand = float(m.group(2).replace(",", "."))
                    continue
    except FileNotFoundError:
        pass

    return {
        "out_path": str(p),
        "run_ms": run_ms,
        "worst_congestion": worst,
        "random_max_congestion": rand,
    }

def resolve_out_path(p, csv_dir):
    # Return absolute path when found so merge keys match
    candidates = [p, os.path.join(csv_dir, p), os.path.join(os.getcwd(), p)]
    for c in candidates:
        if c and os.path.exists(c):
            return os.path.abspath(c)
    return None

def resolve_graph_path(p, csv_dir):
    candidates = []
    candidates.append(p)
    candidates.append(os.path.join(csv_dir, p))
    stem, ext = os.path.splitext(p)
    if ext.lower() != ".lgf":
        candidates.append(stem + ".lgf")
        candidates.append(os.path.join(csv_dir, stem + ".lgf"))
    candidates.append(os.path.join(os.getcwd(), os.path.basename(p)))
    candidates.append(os.path.join(os.getcwd(), os.path.basename(stem) + ".lgf"))
    for sub in ("datasets", "experiments/datasets", "graphs", "logs"):
        candidates.append(os.path.join(sub, os.path.basename(stem) + ".lgf"))
    seen, uniq = set(), []
    for c in candidates:
        if c not in seen:
            uniq.append(c); seen.add(c)
    for c in uniq:
        if c and os.path.exists(c):
            return os.path.abspath(c)
    return None

def count_nodes_lgf(path):
    """
    Supports:
      - DIMACS-like: 'p <type> n m'  -> returns n
      - LEMON .lgf:  '@nodes' section -> counts lines until next '@'
    """
    try:
        with open(path, "r", encoding="utf-8", errors="ignore") as f:
            lines = f.readlines()
    except Exception:
        return None

    # DIMACS 'p' header
    for line in lines:
        s = line.strip()
        if not s or s.startswith(('c','#')):
            continue
        if s.startswith('p '):
            parts = s.split()
            if len(parts) >= 4:
                for cand in (2, -2):
                    try:
                        return int(parts[cand])
                    except Exception:
                        pass
            break

    # LEMON .lgf
    nodes_start = None
    for i, line in enumerate(lines):
        if line.strip().lower().startswith("@nodes"):
            nodes_start = i; break
    if nodes_start is not None:
        i = nodes_start + 1
        while i < len(lines) and not lines[i].strip():  # blank lines
            i += 1
        if i < len(lines) and not lines[i].strip().startswith("@"):
            i += 1  # header row
        cnt = 0
        while i < len(lines):
            s = lines[i].strip()
            if not s:
                i += 1; continue
            if s.startswith("@"):
                break
            cnt += 1; i += 1
        if cnt > 0:
            return cnt
    return None

def grouped_bar(ax, pivot_df, title, ylabel, xlabel="Graph", rotate_xticks=True):
    pivot_df.plot(kind="bar", ax=ax)
    ax.set_title(title, fontsize=10)
    ax.set_ylabel(ylabel, fontsize=9)
    ax.set_xlabel(xlabel, fontsize=2)
    if rotate_xticks:
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
    ax.legend(title="Solver")
    ax.figure.tight_layout()

# ------------------------ Main ------------------------

def main():
    ap = argparse.ArgumentParser(
        description="Plots worst-case & gravity/random TM congestion per graph; |V| vs runtime; plus time summaries."
    )
    ap.add_argument("--csv", required=True, help="Path to experiments CSV (from experiments.py)")
    ap.add_argument("--outdir", default="plots", help="Directory to save plots")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    df = pd.read_csv(args.csv)
    # Add convenience columns if missing
    if "graph_name" not in df.columns:
        df["graph_name"] = df["graph"].apply(normalize_graph_name)
    if "rep" not in df.columns:
        df["rep"] = df["stdout"].apply(extract_rep_from_stdout) if "stdout" in df.columns else None
    if "wall_time_s" in df.columns:
        df["wall_time_s"] = pd.to_numeric(df["wall_time_s"], errors="coerce")

    # If the CSV already has metrics, use them; otherwise parse logs
    need_parse = []
    for col in ["run_ms","worst_congestion","random_max_congestion"]:
        if col not in df.columns or df[col].isna().all():
            need_parse.append(col)

    if need_parse and "stdout" in df.columns:
        csv_dir = os.path.dirname(os.path.abspath(args.csv))
        df["out_abspath"] = df["stdout"].apply(lambda p: resolve_out_path(p, csv_dir))
        out_metrics = []
        for p in df["out_abspath"].dropna().unique():
            out_metrics.append(parse_out_metrics(p))  # accepts str; normalizes + returns abs out_path
        df_out = pd.DataFrame(out_metrics) if out_metrics else pd.DataFrame(
            columns=["out_path","run_ms","worst_congestion","random_max_congestion"]
        )
        dfm = df.merge(df_out, left_on="out_abspath", right_on="out_path", how="left")
        # Prefer CSV values when present; fill from parsed otherwise
        for col in ["run_ms","worst_congestion","random_max_congestion"]:
            if col in df.columns:
                dfm[col] = dfm[col].combine_first(dfm[f"{col}_y"]) if f"{col}_y" in dfm else dfm[col]
        # Drop duplicate suffixes if any
        for col in list(dfm.columns):
            if col.endswith("_x") or col.endswith("_y"):
                base = col[:-2]
                if base in dfm.columns:
                    dfm.drop(columns=[col], inplace=True, errors="ignore")
        # Use dfm going forward
    else:
        dfm = df.copy()

    # Ensure node counts are available
    if "n_nodes" not in dfm.columns or dfm["n_nodes"].isna().all():
        csv_dir = os.path.dirname(os.path.abspath(args.csv))
        dfm["graph_lgf"] = dfm["graph"].apply(lambda p: resolve_graph_path(p, csv_dir))
        dfm["n_nodes"] = dfm["graph_lgf"].apply(lambda p: count_nodes_lgf(p) if p else None)

    # Save merged for debugging
    dfm.to_csv(os.path.join(args.outdir, "merged_with_nodecounts.csv"), index=False)

    # ------------------------ Plots ------------------------

    # Worst-case congestion per graph
    if "worst_congestion" in dfm.columns and dfm["worst_congestion"].notna().any():
        wc = dfm.groupby(["graph_name","solver"], dropna=False)["worst_congestion"].mean().reset_index()
        if not wc.empty:
            pivot = wc.pivot(index="graph_name", columns="solver", values="worst_congestion").sort_index()

            # sorty by cohen
            cohen_col = find_solver_col(pivot.columns, key="cohen")
            if cohen_col is not None:
                order = pivot[cohen_col].sort_values(ascending=True, na_position="last").index
                pivot = pivot.loc[order]
            else:
                pivot = pivot.sort_index()

            fig = plt.figure(); ax = fig.add_subplot(111)
            grouped_bar(ax, pivot, title="Worst-case congestion per graph", ylabel="Worst-case congestion (mean across repeats)")
            plt.savefig(os.path.join(args.outdir, "worst_congestion_per_graph.png"), dpi=180)
            plt.close()

    # Gravity / random traffic matrix max congestion per graph
    if "random_max_congestion" in dfm.columns and dfm["random_max_congestion"].notna().any():
        rc = dfm.groupby(["graph_name","solver"], dropna=False)["random_max_congestion"].mean().reset_index()
        if not rc.empty:
            pivot = rc.pivot(index="graph_name", columns="solver", values="random_max_congestion").sort_index()
            fig = plt.figure(); ax = fig.add_subplot(111)
            grouped_bar(ax, pivot, title="Gravity / Traffic-Matrix Max Congestion per Graph",
                        ylabel="Max congestion on traffic matrix (mean across repeats)")
            plt.savefig(os.path.join(args.outdir, "random_max_congestion_per_graph.png"), dpi=180)
            plt.close()

        # Aggregate worst-case congestion by solver (reusing gravity_max_congestion_by_solver slot)
        agg_r = dfm.groupby("solver", dropna=False)["worst_congestion"].agg(["count","mean","std"]).reset_index()
        if not agg_r.empty:
            fig = plt.figure(); ax = fig.add_subplot(111)
            ax.bar(agg_r["solver"], agg_r["mean"], yerr=agg_r["std"], capsize=4)
            ax.set_title("Worst-case Congestion by Solver")
            ax.set_ylabel("Congestion (mean ± std)")
            ax.set_xlabel("Solver")
            fig.tight_layout()
            plt.savefig(os.path.join(args.outdir, "gravity_max_congestion_by_solver.png"), dpi=180)
            plt.close()


    # Worst vs Gravity scatter
        if "worst_congestion" in dfm.columns and dfm["worst_congestion"].notna().any():
            pair = dfm.groupby(["graph_name","solver"], dropna=False).agg(
                worst=("worst_congestion","mean"),
                rand=("random_max_congestion","mean")
            ).reset_index()
            if not pair.empty:
                fig = plt.figure(); ax = fig.add_subplot(111)
                for solver, sub in pair.groupby("solver"):
                    ax.scatter(sub["worst"], sub["rand"], label=solver)
                ax.set_title("Worst-case vs Gravity / TM Max Congestion")
                ax.set_xlabel("Worst-case congestion (mean)")
                ax.set_ylabel("Gravity / TM max congestion (mean)")
                ax.legend(title="Solver")
                ax.grid(True, which="both", linestyle=":", linewidth=0.5)
                fig.tight_layout()
                plt.savefig(os.path.join(args.outdir, "worst_vs_random_congestion.png"), dpi=180)
                plt.close()

    # |V| vs wall time
    if "n_nodes" in dfm.columns and dfm["n_nodes"].notna().any() and "wall_time_s" in dfm.columns:
        fig = plt.figure(); ax = fig.add_subplot(111)
        for solver, sub in dfm.dropna(subset=["n_nodes","wall_time_s"]).groupby("solver"):
            ax.scatter(sub["n_nodes"], sub["wall_time_s"], label=solver)
        ax.set_title("Number of vertices vs. running time")
        ax.set_xlabel("Number of vertices (|V|)")
        ax.set_ylabel("Wall time (s)")
        ax.legend(title="Solver")
        fig.tight_layout()
        plt.savefig(os.path.join(args.outdir, "nodes_vs_runtime.png"), dpi=180)
        plt.close()

    # Average wall time by solver
    if "wall_time_s" in dfm.columns:
        agg = dfm.groupby("solver", dropna=False)["wall_time_s"].agg(["count","mean","std"]).reset_index()
        if not agg.empty:
            fig = plt.figure(); ax = fig.add_subplot(111)
            ax.bar(agg["solver"], agg["mean"], yerr=agg["std"], capsize=4)
            ax.set_title("Average wall time by solver")
            ax.set_ylabel("Time (s)")
            ax.set_xlabel("Solver")
            fig.tight_layout()
            plt.savefig(os.path.join(args.outdir, "avg_wall_time_by_solver.png"), dpi=180)
            plt.close()

    # Per-run wall times (x-order sorted by number of vertices)
    if "wall_time_s" in dfm.columns and dfm["wall_time_s"].notna().any() and "n_nodes" in dfm.columns:
        runs = dfm.dropna(subset=["wall_time_s", "n_nodes"]).copy()

        # Build a run key to align multiple solvers per instance
        if "rep" in runs.columns and runs["rep"].notna().any():
            runs["run_key"] = runs["graph_name"].astype(str) + "::" + runs["rep"].astype("Int64").astype(str)
        else:
            runs["run_key"] = runs["graph_name"].astype(str)

        # Get mean n_nodes per run_key
        node_order = (
            runs.groupby("run_key")["n_nodes"]
            .mean()
            .sort_values(ascending=True)  # from small to large graphs
            .index
            .tolist()
        )

        xpos = {k: i for i, k in enumerate(node_order)}
        runs["xpos"] = runs["run_key"].map(xpos)

        fig = plt.figure(); ax = fig.add_subplot(111)
        for solver, sub in runs.groupby("solver"):
            sub = sub.sort_values("xpos")
            ax.scatter(sub["xpos"], sub["wall_time_s"], label=solver, marker='o', s=5)

        ax.legend()
        ax.set_title("Per-run wall times (runs ordered by |V|)")
        ax.set_ylabel("Time (s)")
        ax.set_xlabel("Run index (ordered by number of vertices)")
        fig.tight_layout()
        plt.savefig(os.path.join(args.outdir, "per_run_wall_times.png"), dpi=180)
        plt.close()



if __name__ == "__main__":
    main()
