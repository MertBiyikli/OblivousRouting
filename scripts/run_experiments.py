#!/usr/bin/env python3

import argparse
import csv
import json
import subprocess
from datetime import datetime
from pathlib import Path


def safe_name(path_or_name: str) -> str:
    name = Path(path_or_name).stem
    return (
        name.replace(" ", "_")
        .replace("/", "_")
        .replace("\\", "_")
        .replace("(", "")
        .replace(")", "")
    )


def parse_args():
    parser = argparse.ArgumentParser(description="Run oblivious routing experiments.")

    parser.add_argument(
        "--bin",
        required=True,
        help="Path to oblivious_routing binary.",
    )

    parser.add_argument(
        "--graphs",
        nargs="+",
        required=True,
        help="Graph files to run.",
    )

    parser.add_argument(
        "--solvers",
        nargs="+",
        required=True,
        help="Solvers to run, e.g. electrical semi_tree semi_elec.",
    )

    parser.add_argument(
        "--demands",
        default="gravity,bimodal,uniform,gaussian",
        help="Comma-separated demand models.",
    )

    parser.add_argument(
        "--out",
        default=None,
        help="Output directory. Defaults to results/run_<timestamp>.",
    )

    parser.add_argument(
        "--threads",
        type=int,
        default=None,
        help="Optional number of threads.",
    )

    return parser.parse_args()

def as_number(value):
    if value is None:
        return None

    try:
        number = float(value)
        if number.is_integer():
            return int(number)
        return number
    except (TypeError, ValueError):
        return value


def get_nested(data, *keys, default=None):
    current = data
    for key in keys:
        if not isinstance(current, dict):
            return default
        current = current.get(key)
        if current is None:
            return default
    return current

def extract_summary_rows(result_json_path: Path):
    with result_json_path.open("r") as f:
        data = json.load(f)

    rows = []

    # Supports both:
    # 1. current flat format: { "solver": ..., "demand_evaluations": [...] }
    # 2. future nested format: { "solver_results": [...] }
    solver_results = data.get("solver_results")
    if solver_results is None:
        solver_results = [data]

    graph_data = data.get("graph")

    if isinstance(graph_data, dict):
        graph_path = graph_data.get("path")
        graph_name = Path(graph_path).stem if graph_path else None
        nodes = as_number(graph_data.get("nodes") or graph_data.get("num_nodes"))
        edges = as_number(graph_data.get("edges") or graph_data.get("num_edges"))
    else:
        graph_path = graph_data
        graph_name = Path(graph_path).stem if graph_path else None
        nodes = as_number(data.get("nodes"))
        edges = as_number(data.get("edges"))

    for solver_result in solver_results:
        solver_name = (
                solver_result.get("solver")
                or solver_result.get("solver_name")
                or "unknown"
        )

        routing_base = solver_result.get("routing_base")
        status = solver_result.get("status")

        total_runtime = (
                get_nested(solver_result, "runtime", "total_microseconds")
                or solver_result.get("total_runtime_microseconds")
        )

        preprocessing_runtime = (
                get_nested(solver_result, "runtime", "preprocessing_microseconds")
                or solver_result.get("preprocessing_runtime_microseconds")
        )

        solve_runtime = (
                get_nested(solver_result, "runtime", "solve_microseconds")
                or solver_result.get("solve_runtime_microseconds")
        )

        candidate_paths = (
                get_nested(solver_result, "routing_scheme", "candidate_paths")
                or solver_result.get("candidate_paths")
        )

        avg_paths = (
                get_nested(solver_result, "routing_scheme", "average_paths_per_pair")
                or solver_result.get("average_paths_per_pair")
        )

        top_level_oblivious_ratio = solver_result.get("oblivious_ratio")

        mwu = solver_result.get("mwu_metrics") or {}

        for eval_result in solver_result.get("demand_evaluations", []):
            rows.append({
                "graph": graph_name,
                "graph_path": graph_path,
                "nodes": nodes,
                "edges": edges,

                "solver": solver_name,
                "routing_base": routing_base,
                "status": status,

                "demand_model": eval_result.get("demand_model"),
                "congestion": as_number(eval_result.get("congestion")),
                "offline_opt": as_number(eval_result.get("offline_opt")),
                "oblivious_ratio": as_number(
                    eval_result.get("oblivious_ratio", top_level_oblivious_ratio)
                ),

                "total_runtime_microseconds": as_number(total_runtime),
                "preprocessing_runtime_microseconds": as_number(preprocessing_runtime),
                "solve_runtime_microseconds": as_number(solve_runtime),
                "evaluation_runtime_microseconds": as_number(
                    eval_result.get("runtime_microseconds")
                ),

                "candidate_paths": as_number(candidate_paths),
                "average_paths_per_pair": as_number(avg_paths),

                "mwu_iteration_count": as_number(mwu.get("iteration_count")),
                "mwu_solve_time_microseconds": as_number(mwu.get("solve_time_microseconds")),
                "mwu_transformation_time_microseconds": as_number(
                    mwu.get("transformation_time_microseconds")
                ),
                "mwu_load_computation_time_microseconds": as_number(
                    mwu.get("load_computation_time_microseconds")
                ),
                "mwu_weight_update_time_microseconds": as_number(
                    mwu.get("weight_update_time_microseconds")
                ),
                "mwu_average_oracle_time_microseconds": as_number(
                    mwu.get("average_oracle_time_microseconds")
                ),

                "json_file": str(result_json_path),
            })

    return rows


def main():
    args = parse_args()

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_dir = Path(args.out) if args.out else Path("results") / f"run_{timestamp}"
    out_dir.mkdir(parents=True, exist_ok=True)

    all_rows = []
    failures = []

    for graph in args.graphs:
        for solver in args.solvers:
            graph_name = safe_name(graph)
            solver_name = safe_name(solver)

            json_path = out_dir / f"{graph_name}__{solver_name}.json"
            stderr_path = out_dir / f"{graph_name}__{solver_name}.stderr.txt"

            cmd = [
                args.bin,
                solver,
                graph,
                args.demands,
                str(json_path),
            ]

            if args.threads is not None:
                cmd += ["--threads", str(args.threads)]

            print(f"[RUN] {graph_name} | {solver}")

            completed = subprocess.run(
                cmd,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )

            if completed.returncode != 0:
                stderr_path.write_text(completed.stderr)
                failures.append({
                    "graph": graph,
                    "solver": solver,
                    "returncode": completed.returncode,
                    "stderr_file": str(stderr_path),
                })
                print(f"[FAIL] {graph_name} | {solver}")
                continue

            try:
                rows = extract_summary_rows(json_path)
                all_rows.extend(rows)
                print(f"[OK]   {graph_name} | {solver}")
            except Exception as e:
                stderr_path.write_text(
                    f"JSON parsing/extraction failed:\n{e}\n\n"
                    f"stdout:\n{completed.stdout}\n\n"
                    f"stderr:\n{completed.stderr}\n"
                )
                failures.append({
                    "graph": graph,
                    "solver": solver,
                    "returncode": "json_error",
                    "stderr_file": str(stderr_path),
                })
                print(f"[BAD JSON] {graph_name} | {solver}")

    summary_path = out_dir / "summary.csv"

    fieldnames = [
        "graph",
        "graph_path",
        "nodes",
        "edges",

        "solver",
        "routing_base",
        "status",

        "demand_model",
        "congestion",
        "offline_opt",
        "oblivious_ratio",

        "total_runtime_microseconds",
        "preprocessing_runtime_microseconds",
        "solve_runtime_microseconds",
        "evaluation_runtime_microseconds",

        "candidate_paths",
        "average_paths_per_pair",

        "mwu_iteration_count",
        "mwu_solve_time_microseconds",
        "mwu_transformation_time_microseconds",
        "mwu_load_computation_time_microseconds",
        "mwu_weight_update_time_microseconds",
        "mwu_average_oracle_time_microseconds",

        "json_file",
    ]

    with summary_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(all_rows)

    if failures:
        failure_path = out_dir / "failures.json"
        failure_path.write_text(json.dumps(failures, indent=2))
        print()
        print(f"Finished with {len(failures)} failed runs.")
        print(f"Failures written to: {failure_path}")

    print()
    print(f"Results directory: {out_dir}")
    print(f"Summary CSV:       {summary_path}")

    if failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()