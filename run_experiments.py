#!/usr/bin/env python3
import subprocess
import re
import csv
import json
import argparse
import glob
from pathlib import Path
from typing import List, Dict


# ======================
# REGEX DEFINITIONS
# ======================

RE_GRAPH = re.compile(r"Graph loaded:\s+(\d+)\s+nodes,\s+(\d+)\s+edges")
RE_SOLVER_HEADER = re.compile(r"=== Running solver: (.+?) ===")
RE_TOTAL_TIME = re.compile(r"Total running time:\s+(\d+)\s+ms")
RE_SOLVE_TIME = re.compile(r"Solve time:\s+(\d+)\s+ms")
RE_TRANSFORM_TIME = re.compile(r"Transformation time:\s+(\d+)\s+ms")
RE_MWU_ITERS = re.compile(r"MWU iterations:\s+(\d+)")
RE_AVG_ORACLE = re.compile(r"Average oracle time:\s+([\d.]+)\s+ms")
RE_RATIO = re.compile(
    r"Ratio off the optimal offline solution .*?:\s+([\d.]+)%\s+\(([\d.]+)\s+/\s+([\d.]+)\)"
)


# ======================
# PARSING
# ======================

def parse_output(stdout: str, dataset: str, demand_model: str) -> List[Dict]:
    results = []
    current = None

    num_nodes = None
    num_edges = None

    for line in stdout.splitlines():

        # ---- Graph metadata (once per run) ----
        if m := RE_GRAPH.search(line):
            num_nodes = int(m.group(1))
            num_edges = int(m.group(2))

        # ---- Solver-specific metrics----
        if m := RE_SOLVER_HEADER.search(line):
            if current:
                results.append(current)
            current = {
                "dataset": Path(dataset).name,
                "demand_model": demand_model,
                "solver": m.group(1),
                "num_nodes": num_nodes,
                "num_edges": num_edges,
            }

        elif current is not None:
            if m := RE_TOTAL_TIME.search(line):
                current["total_time_ms"] = int(m.group(1))
            elif m := RE_SOLVE_TIME.search(line):
                current["solve_time_ms"] = int(m.group(1))
            elif m := RE_TRANSFORM_TIME.search(line):
                current["transform_time_ms"] = int(m.group(1))
            elif m := RE_MWU_ITERS.search(line):
                current["mwu_iterations"] = int(m.group(1))
            elif m := RE_AVG_ORACLE.search(line):
                current["avg_oracle_time_ms"] = float(m.group(1))
            elif m := RE_RATIO.search(line):
                current["ratio_percent"] = float(m.group(1))
                current["offline_opt_value"] = float(m.group(2))
                current["achieved_congestion"] = float(m.group(3))

    if current:
        results.append(current)

    return results


# ======================
# RUN BINARY
# ======================

def run_binary(binary: str, solvers: str, dataset: str, demand: str) -> str:
    cmd = [binary, solvers, dataset, demand]

    print("▶ Running:", " ".join(cmd))
    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        check=True,
    )
    return result.stdout


# ======================
# MAIN
# ======================

def main():
    parser = argparse.ArgumentParser(
        description="Run Oblivious Routing experiments"
    )

    parser.add_argument(
        "--solvers",
        required=True,
        help="Comma-separated solver IDs, e.g. 0,1,2"
    )

    parser.add_argument(
        "--datasets",
        required=True,
        nargs="+",
        help="List of .lgf datasets"
    )

    parser.add_argument(
        "--demands",
        nargs="+",
        default=["gravity"],
        help="Demand models to evaluate (e.g. gravity uniform poisson)"
    )


    parser.add_argument(
        "--binary",
        default="./cmake-build-release/oblivious_routing",
        help="Path to oblivious_routing binary"
    )

    parser.add_argument(
        "--out",
        default="results",
        help="Output directory"
    )

    args = parser.parse_args()

    out_dir = Path(args.out)
    out_dir.mkdir(exist_ok=True)

    all_results = []

    expanded_datasets = []

    for pattern in args.datasets:
        matches = glob.glob(pattern)
        if not matches:
            raise RuntimeError(f"No datasets matched pattern: {pattern}")
        expanded_datasets.extend(matches)

    expanded_datasets = sorted(set(expanded_datasets))

    for dataset in expanded_datasets:
        for demand in args.demands:
            stdout = run_binary(
                binary=args.binary,
                solvers=args.solvers,
                dataset=dataset,
                demand=demand,
            )

            results = parse_output(
                stdout=stdout,
                dataset=dataset,
                demand_model=demand,
            )

            all_results.extend(results)


    # ======================
    # SAVE RESULTS
    # ======================
    csv_file = out_dir / "results.csv"
    json_file = out_dir / "results.json"

    if all_results:
        with csv_file.open("w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=all_results[0].keys())
            writer.writeheader()
            writer.writerows(all_results)

        with json_file.open("w") as f:
            json.dump(all_results, f, indent=2)

    print(f"\n✔ Saved {len(all_results)} records")
    print(f"  CSV : {csv_file}")
    print(f"  JSON: {json_file}")


if __name__ == "__main__":
    main()
