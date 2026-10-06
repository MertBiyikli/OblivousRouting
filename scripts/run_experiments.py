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
    parser = argparse.ArgumentParser(
        description="Run E-Routing experiments."
    )

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
        help=(
            "Solvers to run, e.g. "
            "electrical raecke_ckr expander."
        ),
    )

    parser.add_argument(
        "--demands",
        default="gravity,bimodal,uniform,gaussian",
        help="Comma-separated demand models.",
    )

    parser.add_argument(
        "--out",
        default=None,
        help=(
            "Output directory. Defaults to "
            "results/run_<timestamp>."
        ),
    )

    parser.add_argument(
        "--threads",
        type=int,
        default=None,
        help="Optional number of worker threads.",
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


def extract_summary_rows(result_json_path: Path):
    """
    Convert one schema-1.0 E-Routing result JSON into
    normalized CSV rows.

    One CSV row represents:

        graph x solver x demand model
    """

    with result_json_path.open() as f:
        data = json.load(f)

    schema_version = data.get("schema_version")

    if schema_version != "1.0":
        raise ValueError(
            f"Unsupported result schema version: "
            f"{schema_version!r}"
        )

    # ------------------------------------------------------------------
    # Graph metadata
    # ------------------------------------------------------------------

    graph = data.get("graph", {})

    if not isinstance(graph, dict):
        raise ValueError(
            "Invalid schema: 'graph' must be an object."
        )

    graph_name = graph.get("name")
    graph_path = graph.get("path")

    nodes = as_number(
        graph.get("nodes")
    )

    edges = as_number(
        graph.get("edges")
    )

    # ------------------------------------------------------------------
    # Experiment configuration
    # ------------------------------------------------------------------

    configuration = data.get(
        "configuration",
        {}
    )

    if not isinstance(configuration, dict):
        raise ValueError(
            "Invalid schema: "
            "'configuration' must be an object."
        )

    seed = as_number(
        configuration.get("seed")
    )

    threads = as_number(
        configuration.get("threads")
    )

    # ------------------------------------------------------------------
    # Solver results
    # ------------------------------------------------------------------

    solver_results = data.get(
        "solver_results"
    )

    if not isinstance(solver_results, list):
        raise ValueError(
            "Invalid schema: "
            "'solver_results' must be an array."
        )

    rows = []

    for solver_result in solver_results:

        if not isinstance(solver_result, dict):
            raise ValueError(
                "Invalid schema: "
                "solver result must be an object."
            )

        solver_name = solver_result.get(
            "solver"
        )

        solver_type = solver_result.get(
            "solver_type"
        )

        status = solver_result.get(
            "status"
        )

        routing_base = solver_result.get(
            "routing_base"
        )

        # --------------------------------------------------------------
        # Common solver metrics
        # --------------------------------------------------------------

        runtime = solver_result.get(
            "runtime",
            {}
        )

        quality = solver_result.get(
            "quality",
            {}
        )

        routing_scheme = solver_result.get(
            "routing_scheme",
            {}
        )

        algorithm_metrics = solver_result.get(
            "algorithm_metrics",
            {}
        )

        if not isinstance(runtime, dict):
            runtime = {}

        if not isinstance(quality, dict):
            quality = {}

        if not isinstance(routing_scheme, dict):
            routing_scheme = {}

        if not isinstance(algorithm_metrics, dict):
            algorithm_metrics = {}

        total_runtime = as_number(
            runtime.get(
                "total_microseconds"
            )
        )

        preprocessing_runtime = as_number(
            runtime.get(
                "preprocessing_microseconds"
            )
        )

        solve_runtime = as_number(
            runtime.get(
                "solve_microseconds"
            )
        )

        oblivious_ratio = as_number(
            quality.get(
                "oblivious_ratio"
            )
        )

        candidate_paths = as_number(
            routing_scheme.get(
                "candidate_paths"
            )
        )

        average_paths_per_pair = as_number(
            routing_scheme.get(
                "average_paths_per_pair"
            )
        )

        # --------------------------------------------------------------
        # Solver-specific metrics
        # --------------------------------------------------------------

        mwu = algorithm_metrics.get(
            "mwu",
            {}
        )

        expander = algorithm_metrics.get(
            "expander",
            {}
        )

        if not isinstance(mwu, dict):
            mwu = {}

        if not isinstance(expander, dict):
            expander = {}

        # --------------------------------------------------------------
        # Demand evaluations
        # --------------------------------------------------------------

        demand_evaluations = solver_result.get(
            "demand_evaluations",
            []
        )

        if not isinstance(demand_evaluations, list):
            raise ValueError(
                "Invalid schema: "
                "'demand_evaluations' must be an array."
            )

        for evaluation in demand_evaluations:

            if not isinstance(evaluation, dict):
                continue

            rows.append({

                # Schema
                "schema_version":
                    schema_version,

                # Graph
                "graph":
                    graph_name,

                "graph_path":
                    graph_path,

                "nodes":
                    nodes,

                "edges":
                    edges,

                # Configuration
                "seed":
                    seed,

                "threads":
                    threads,

                # Solver
                "solver":
                    solver_name,

                "solver_type":
                    solver_type,

                "routing_base":
                    routing_base,

                "status":
                    status,

                # Demand result
                "demand_model":
                    evaluation.get(
                        "demand_model"
                    ),

                "congestion":
                    as_number(
                        evaluation.get(
                            "congestion"
                        )
                    ),

                "evaluation_runtime_microseconds":
                    as_number(
                        evaluation.get(
                            "runtime_microseconds"
                        )
                    ),

                # Common quality
                "oblivious_ratio":
                    oblivious_ratio,

                # Runtime
                "total_runtime_microseconds":
                    total_runtime,

                "preprocessing_runtime_microseconds":
                    preprocessing_runtime,

                "solve_runtime_microseconds":
                    solve_runtime,

                # Routing scheme
                "candidate_paths":
                    candidate_paths,

                "average_paths_per_pair":
                    average_paths_per_pair,

                # ------------------------------------------------------
                # MWU metrics
                # ------------------------------------------------------

                "mwu_iteration_count":
                    as_number(
                        mwu.get(
                            "iteration_count"
                        )
                    ),

                "mwu_solve_time_microseconds":
                    as_number(
                        mwu.get(
                            "solve_time_microseconds"
                        )
                    ),

                "mwu_transformation_time_microseconds":
                    as_number(
                        mwu.get(
                            "transformation_time_microseconds"
                        )
                    ),

                "mwu_load_computation_time_microseconds":
                    as_number(
                        mwu.get(
                            "load_computation_time_microseconds"
                        )
                    ),

                "mwu_weight_update_time_microseconds":
                    as_number(
                        mwu.get(
                            "weight_update_time_microseconds"
                        )
                    ),

                "mwu_average_oracle_time_microseconds":
                    as_number(
                        mwu.get(
                            "average_oracle_time_microseconds"
                        )
                    ),

                "mwu_oracle_calls":
                    as_number(
                        mwu.get(
                            "oracle_calls"
                        )
                    ),

                # ------------------------------------------------------
                # Expander metrics
                # ------------------------------------------------------

                "expander_hierarchy_runtime_microseconds":
                    as_number(
                        expander.get(
                            "hierarchy_runtime_microseconds"
                        )
                    ),

                "expander_tree_runtime_microseconds":
                    as_number(
                        expander.get(
                            "tree_runtime_microseconds"
                        )
                    ),

                "expander_basis_flow_runtime_microseconds":
                    as_number(
                        expander.get(
                            "basis_flow_runtime_microseconds"
                        )
                    ),

                "expander_hierarchy_levels":
                    as_number(
                        expander.get(
                            "hierarchy_levels"
                        )
                    ),

                "expander_hierarchy_clusters":
                    as_number(
                        expander.get(
                            "hierarchy_clusters"
                        )
                    ),

                "expander_tree_nodes":
                    as_number(
                        expander.get(
                            "tree_nodes"
                        )
                    ),

                "expander_tree_edges":
                    as_number(
                        expander.get(
                            "tree_edges"
                        )
                    ),

                "expander_tree_depth":
                    as_number(
                        expander.get(
                            "tree_depth"
                        )
                    ),

                "expander_basis_flows":
                    as_number(
                        expander.get(
                            "basis_flows"
                        )
                    ),

                "expander_total_electrical_solves":
                    as_number(
                        expander.get(
                            "total_electrical_solves"
                        )
                    ),

                "expander_average_electrical_solves":
                    as_number(
                        expander.get(
                            "average_electrical_solves"
                        )
                    ),

                "expander_max_basis_embedding_congestion":
                    as_number(
                        expander.get(
                            "max_basis_embedding_congestion"
                        )
                    ),

                "expander_max_conservation_error":
                    as_number(
                        expander.get(
                            "max_conservation_error"
                        )
                    ),

                # Source result file
                "json_file":
                    str(result_json_path),
            })

    return rows


def main():
    args = parse_args()

    timestamp = datetime.now().strftime(
        "%Y%m%d_%H%M%S"
    )

    out_dir = (
        Path(args.out)
        if args.out
        else Path("results") / f"run_{timestamp}"
    )

    out_dir.mkdir(
        parents=True,
        exist_ok=True
    )

    all_rows = []
    failures = []

    # ------------------------------------------------------------------
    # The C++ CLI accepts comma-separated solvers:
    #
    #   --solver electrical,raecke_ckr,expander
    #
    # Therefore we run the binary once per graph instead of once for
    # every graph/solver combination.
    # ------------------------------------------------------------------

    solver_argument = ",".join(
        args.solvers
    )

    for graph in args.graphs:

        graph_name = safe_name(
            graph
        )

        json_path = (
                out_dir /
                f"{graph_name}.json"
        )

        stderr_path = (
                out_dir /
                f"{graph_name}.stderr.txt"
        )

        cmd = [
            args.bin,
            "solve",

            "--solver",
            solver_argument,

            "--graph",
            graph,

            "--demand",
            args.demands,

            "--output",
            str(json_path),

            "--output-format",
            "json",
        ]

        if args.threads is not None:
            cmd += [
                "--threads",
                str(args.threads),
            ]

        print(
            f"[RUN] {graph_name}"
        )

        print(
            f"      solvers: "
            f"{solver_argument}"
        )

        print(
            f"      demands: "
            f"{args.demands}"
        )

        completed = subprocess.run(
            cmd,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )

        # --------------------------------------------------------------
        # C++ execution failed
        # --------------------------------------------------------------

        if completed.returncode != 0:

            stderr_path.write_text(
                completed.stderr
            )

            failures.append({
                "graph":
                    graph,

                "solvers":
                    args.solvers,

                "returncode":
                    completed.returncode,

                "stderr_file":
                    str(stderr_path),
            })

            print(
                f"[FAIL] {graph_name}"
            )

            continue

        # --------------------------------------------------------------
        # Parse standardized result
        # --------------------------------------------------------------

        try:

            rows = extract_summary_rows(
                json_path
            )

            all_rows.extend(
                rows
            )

            print(
                f"[OK]   {graph_name} "
                f"| {len(rows)} summary rows"
            )

        except Exception as exc:

            stderr_path.write_text(
                "JSON parsing/extraction failed:\n"
                f"{exc}\n\n"
                f"stdout:\n"
                f"{completed.stdout}\n\n"
                f"stderr:\n"
                f"{completed.stderr}\n"
            )

            failures.append({
                "graph":
                    graph,

                "solvers":
                    args.solvers,

                "returncode":
                    "json_error",

                "stderr_file":
                    str(stderr_path),
            })

            print(
                f"[BAD JSON] {graph_name}"
            )

    # ------------------------------------------------------------------
    # Summary CSV
    # ------------------------------------------------------------------

    summary_path = (
            out_dir /
            "summary.csv"
    )

    fieldnames = [

        "schema_version",

        # Graph
        "graph",
        "graph_path",
        "nodes",
        "edges",

        # Configuration
        "seed",
        "threads",

        # Solver
        "solver",
        "solver_type",
        "routing_base",
        "status",

        # Demand evaluation
        "demand_model",
        "congestion",
        "oblivious_ratio",

        # Runtime
        "total_runtime_microseconds",
        "preprocessing_runtime_microseconds",
        "solve_runtime_microseconds",
        "evaluation_runtime_microseconds",

        # Routing scheme
        "candidate_paths",
        "average_paths_per_pair",

        # MWU
        "mwu_iteration_count",
        "mwu_solve_time_microseconds",
        "mwu_transformation_time_microseconds",
        "mwu_load_computation_time_microseconds",
        "mwu_weight_update_time_microseconds",
        "mwu_average_oracle_time_microseconds",
        "mwu_oracle_calls",

        # Expander hierarchy
        "expander_hierarchy_runtime_microseconds",
        "expander_tree_runtime_microseconds",
        "expander_basis_flow_runtime_microseconds",
        "expander_hierarchy_levels",
        "expander_hierarchy_clusters",
        "expander_tree_nodes",
        "expander_tree_edges",
        "expander_tree_depth",
        "expander_basis_flows",
        "expander_total_electrical_solves",
        "expander_average_electrical_solves",
        "expander_max_basis_embedding_congestion",
        "expander_max_conservation_error",

        # Origin
        "json_file",
    ]

    with summary_path.open(
            "w",
            newline=""
    ) as f:

        writer = csv.DictWriter(
            f,
            fieldnames=fieldnames
        )

        writer.writeheader()

        writer.writerows(
            all_rows
        )

    # ------------------------------------------------------------------
    # Failures
    # ------------------------------------------------------------------

    if failures:

        failure_path = (
                out_dir /
                "failures.json"
        )

        failure_path.write_text(
            json.dumps(
                failures,
                indent=2
            )
        )

        print()

        print(
            f"Finished with "
            f"{len(failures)} failed graph runs."
        )

        print(
            f"Failures written to: "
            f"{failure_path}"
        )

    # ------------------------------------------------------------------
    # Final report
    # ------------------------------------------------------------------

    print()

    print(
        f"Results directory: "
        f"{out_dir}"
    )

    print(
        f"Summary CSV:       "
        f"{summary_path}"
    )

    print(
        f"Summary rows:      "
        f"{len(all_rows)}"
    )

    if failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()