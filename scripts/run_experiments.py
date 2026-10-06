#!/usr/bin/env python3

import argparse
import csv
import hashlib
import json
import subprocess
from collections import defaultdict
from datetime import datetime
from pathlib import Path

from generate_report import generate_report_from_summary


GRAPH_SUFFIXES = {
    ".lgf",
}


# ============================================================================
# Helpers
# ============================================================================


def safe_name(value):
    name = str(value)

    for char in (
            " ",
            "/",
            "\\",
            "(",
            ")",
            ":",
    ):
        name = name.replace(
            char,
            "_",
        )

    return (
            name.strip("_")
            or "graph"
    )


def as_number(value):
    if (
            value is None
            or value == ""
    ):
        return None

    try:
        number = float(
            value
        )

        if number.is_integer():
            return int(
                number
            )

        return number

    except (
            TypeError,
            ValueError,
    ):
        return value


# ============================================================================
# CLI
# ============================================================================


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Run E-Routing experiments, generate summary.csv "
            "and automatically generate benchmark reports."
        )
    )

    parser.add_argument(
        "--bin",
        required=True,
        help=(
            "Path to oblivious_routing binary."
        ),
    )

    parser.add_argument(
        "--graphs",
        nargs="+",
        required=True,
        help=(
            "Graph files and/or directories. "
            "Directories are searched recursively "
            "for .lgf files."
        ),
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
        default=(
            "gravity,bimodal,"
            "uniform,gaussian"
        ),
        help=(
            "Comma-separated demand models."
        ),
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
        help=(
            "Optional number of worker threads."
        ),
    )

    parser.add_argument(
        "--report-title",
        default=(
            "E-Routing Benchmark Report"
        ),
        help=(
            "Title used for the generated report."
        ),
    )

    return parser.parse_args()


# ============================================================================
# Graph discovery
# ============================================================================


def discover_graphs(inputs):
    discovered = []

    for raw in inputs:
        path = Path(
            raw
        ).expanduser()

        if not path.exists():
            raise ValueError(
                f"Graph input does not exist: "
                f"{path}"
            )

        if path.is_file():

            if (
                    path.suffix.lower()
                    not in GRAPH_SUFFIXES
            ):
                raise ValueError(
                    f"Unsupported graph file: "
                    f"{path}"
                )

            discovered.append(
                path
            )

            continue

        if path.is_dir():

            matches = sorted(
                p
                for p in path.rglob("*")
                if (
                        p.is_file()
                        and
                        p.suffix.lower()
                        in GRAPH_SUFFIXES
                )
            )

            discovered.extend(
                matches
            )

            continue

        raise ValueError(
            f"Unsupported graph input: "
            f"{path}"
        )

    unique = []
    seen = set()

    for path in discovered:
        resolved = (
            path.resolve()
        )

        if resolved in seen:
            continue

        seen.add(
            resolved
        )

        unique.append(
            path
        )

    if not unique:
        raise ValueError(
            "No .lgf graph files were discovered."
        )

    return unique


def make_output_stems(
        graph_paths,
):
    by_stem = defaultdict(
        list
    )

    for path in graph_paths:
        by_stem[
            path.stem
        ].append(path)

    result = {}

    for (
            stem,
            paths,
    ) in by_stem.items():

        if len(paths) == 1:

            result[
                paths[0]
            ] = safe_name(
                stem
            )

            continue

        for path in paths:

            digest = hashlib.sha1(
                str(
                    path.resolve()
                ).encode(
                    "utf-8"
                )
            ).hexdigest()[:8]

            result[path] = (
                f"{safe_name(stem)}"
                f"__{digest}"
            )

    return result


# ============================================================================
# JSON -> summary.csv
# ============================================================================


def extract_summary_rows(
        result_json_path,
):
    with Path(
            result_json_path
    ).open(
        encoding="utf-8"
    ) as f:

        data = json.load(
            f
        )

    schema_version = (
        data.get(
            "schema_version"
        )
    )

    if (
            schema_version
            != "1.0"
    ):
        raise ValueError(
            "Unsupported result schema "
            f"version: {schema_version!r}"
        )

    graph = data.get(
        "graph",
        {},
    )

    configuration = data.get(
        "configuration",
        {},
    )

    solver_results = data.get(
        "solver_results"
    )

    if not isinstance(
            graph,
            dict,
    ):
        raise ValueError(
            "'graph' must be an object."
        )

    if not isinstance(
            configuration,
            dict,
    ):
        raise ValueError(
            "'configuration' must be an object."
        )

    if not isinstance(
            solver_results,
            list,
    ):
        raise ValueError(
            "'solver_results' must be an array."
        )

    rows = []

    for solver_result in solver_results:

        runtime = (
            solver_result.get(
                "runtime",
                {},
            )
        )

        quality = (
            solver_result.get(
                "quality",
                {},
            )
        )

        routing_scheme = (
            solver_result.get(
                "routing_scheme",
                {},
            )
        )

        metrics = (
            solver_result.get(
                "algorithm_metrics",
                {},
            )
        )

        if not isinstance(
                runtime,
                dict,
        ):
            runtime = {}

        if not isinstance(
                quality,
                dict,
        ):
            quality = {}

        if not isinstance(
                routing_scheme,
                dict,
        ):
            routing_scheme = {}

        if not isinstance(
                metrics,
                dict,
        ):
            metrics = {}

        mwu = metrics.get(
            "mwu",
            {},
        )

        expander = metrics.get(
            "expander",
            {},
        )

        if not isinstance(
                mwu,
                dict,
        ):
            mwu = {}

        if not isinstance(
                expander,
                dict,
        ):
            expander = {}

        common = {
            "schema_version":
                schema_version,

            "graph":
                graph.get(
                    "name"
                ),

            "graph_path":
                graph.get(
                    "path"
                ),

            "nodes":
                as_number(
                    graph.get(
                        "nodes"
                    )
                ),

            "edges":
                as_number(
                    graph.get(
                        "edges"
                    )
                ),

            "seed":
                as_number(
                    configuration.get(
                        "seed"
                    )
                ),

            "threads":
                as_number(
                    configuration.get(
                        "threads"
                    )
                ),

            "solver":
                solver_result.get(
                    "solver"
                ),

            "solver_type":
                solver_result.get(
                    "solver_type"
                ),

            "routing_base":
                solver_result.get(
                    "routing_base"
                ),

            "status":
                solver_result.get(
                    "status"
                ),

            "oblivious_ratio":
                as_number(
                    quality.get(
                        "oblivious_ratio"
                    )
                ),

            "total_runtime_microseconds":
                as_number(
                    runtime.get(
                        "total_microseconds"
                    )
                ),

            "preprocessing_runtime_microseconds":
                as_number(
                    runtime.get(
                        "preprocessing_microseconds"
                    )
                ),

            "solve_runtime_microseconds":
                as_number(
                    runtime.get(
                        "solve_microseconds"
                    )
                ),

            "candidate_paths":
                as_number(
                    routing_scheme.get(
                        "candidate_paths"
                    )
                ),

            "average_paths_per_pair":
                as_number(
                    routing_scheme.get(
                        "average_paths_per_pair"
                    )
                ),

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

            "json_file":
                str(
                    result_json_path
                ),
        }

        demand_evaluations = (
            solver_result.get(
                "demand_evaluations",
                [],
            )
        )

        if not isinstance(
                demand_evaluations,
                list,
        ):
            raise ValueError(
                "'demand_evaluations' "
                "must be an array."
            )

        for evaluation in demand_evaluations:

            if not isinstance(
                    evaluation,
                    dict,
            ):
                continue

            row = dict(
                common
            )

            row.update({
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
            })

            rows.append(
                row
            )

    return rows


SUMMARY_FIELDS = [
    "schema_version",

    "graph",
    "graph_path",
    "nodes",
    "edges",

    "seed",
    "threads",

    "solver",
    "solver_type",
    "routing_base",
    "status",

    "demand_model",
    "congestion",
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
    "mwu_oracle_calls",

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

    "json_file",
]


def write_summary(
        summary_path,
        rows,
):
    with Path(
            summary_path
    ).open(
        "w",
        newline="",
        encoding="utf-8",
    ) as f:

        writer = csv.DictWriter(
            f,
            fieldnames=SUMMARY_FIELDS,
        )

        writer.writeheader()

        writer.writerows(
            rows
        )


# ============================================================================
# Main
# ============================================================================


def main():
    args = parse_args()

    if (
            args.threads is not None
            and args.threads <= 0
    ):
        raise SystemExit(
            "--threads must be positive."
        )

    binary = Path(
        args.bin
    ).expanduser()

    if not binary.exists():
        raise SystemExit(
            f"Binary does not exist: "
            f"{binary}"
        )

    try:
        graph_paths = (
            discover_graphs(
                args.graphs
            )
        )

    except ValueError as exc:
        raise SystemExit(
            str(exc)
        ) from exc

    timestamp = (
        datetime.now().strftime(
            "%Y%m%d_%H%M%S"
        )
    )

    out_dir = (
        Path(
            args.out
        ).expanduser()
        if args.out
        else
        Path(
            "results"
        )
        / f"run_{timestamp}"
    )

    out_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    output_stems = (
        make_output_stems(
            graph_paths
        )
    )

    solver_argument = ",".join(
        args.solvers
    )

    all_rows = []
    failures = []

    print()
    print(
        f"Discovered "
        f"{len(graph_paths)} graph(s)."
    )

    print(
        f"Results directory: "
        f"{out_dir}"
    )

    print()

    # ------------------------------------------------------------------
    # Run experiments
    # ------------------------------------------------------------------

    for (
            index,
            graph_path,
    ) in enumerate(
        graph_paths,
        start=1,
    ):

        stem = (
            output_stems[
                graph_path
            ]
        )

        json_path = (
                out_dir
                / f"{stem}.json"
        )

        stderr_path = (
                out_dir
                / f"{stem}.stderr.txt"
        )

        cmd = [
            str(binary),

            "solve",

            "--solver",
            solver_argument,

            "--graph",
            str(graph_path),

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
                str(
                    args.threads
                ),
            ]

        print(
            f"[RUN {index}/"
            f"{len(graph_paths)}] "
            f"{graph_path}"
        )

        completed = (
            subprocess.run(
                cmd,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
        )

        if (
                completed.returncode
                != 0
        ):

            stderr_path.write_text(
                completed.stderr,
                encoding="utf-8",
            )

            failures.append({
                "graph":
                    str(
                        graph_path
                    ),

                "solvers":
                    args.solvers,

                "returncode":
                    completed.returncode,

                "stderr_file":
                    str(
                        stderr_path
                    ),
            })

            print(
                f"[FAIL] "
                f"{graph_path}"
            )

            continue

        try:

            rows = (
                extract_summary_rows(
                    json_path
                )
            )

            all_rows.extend(
                rows
            )

            print(
                f"[OK]   "
                f"{graph_path} "
                f"| {len(rows)} rows"
            )

        except Exception as exc:

            stderr_path.write_text(
                (
                    "JSON parsing/extraction failed:\n"
                    f"{exc}\n\n"
                    f"stdout:\n"
                    f"{completed.stdout}\n\n"
                    f"stderr:\n"
                    f"{completed.stderr}\n"
                ),
                encoding="utf-8",
            )

            failures.append({
                "graph":
                    str(
                        graph_path
                    ),

                "solvers":
                    args.solvers,

                "returncode":
                    "json_error",

                "stderr_file":
                    str(
                        stderr_path
                    ),
            })

            print(
                f"[BAD JSON] "
                f"{graph_path}"
            )

    # ------------------------------------------------------------------
    # summary.csv
    # ------------------------------------------------------------------

    summary_path = (
            out_dir
            / "summary.csv"
    )

    write_summary(
        summary_path,
        all_rows,
    )

    # ------------------------------------------------------------------
    # Report
    # ------------------------------------------------------------------

    report_path = (
            out_dir
            / "report.md"
    )

    report_error = None
    report_outputs = None

    if all_rows:

        try:

            report_outputs = (
                generate_report_from_summary(
                    summary_path=summary_path,
                    report_path=report_path,
                    title=args.report_title,
                )
            )

            print()
            print(
                "[REPORT] Generated:"
            )

            print(
                f"  Markdown: "
                f"{report_outputs['markdown']}"
            )

            if (
                    report_outputs[
                        "html"
                    ]
                    is not None
            ):

                print(
                    f"  HTML:     "
                    f"{report_outputs['html']}"
                )

            else:

                print(
                    "  HTML:     skipped"
                )

            if (
                    report_outputs[
                        "pdf"
                    ]
                    is not None
            ):

                print(
                    f"  PDF:      "
                    f"{report_outputs['pdf']}"
                )

            else:

                print(
                    "  PDF:      skipped"
                )

        except Exception as exc:

            report_error = str(
                exc
            )

            print()
            print(
                f"[REPORT ERROR] "
                f"{report_error}"
            )

    else:

        report_error = (
            "No successful experiment "
            "rows were produced."
        )

    # ------------------------------------------------------------------
    # failures.json
    # ------------------------------------------------------------------

    if failures:

        failure_path = (
                out_dir
                / "failures.json"
        )

        failure_path.write_text(
            json.dumps(
                failures,
                indent=2,
            ),
            encoding="utf-8",
        )

        print()
        print(
            f"Finished with "
            f"{len(failures)} "
            f"failed graph run(s)."
        )

        print(
            f"Failures: "
            f"{failure_path}"
        )

    # ------------------------------------------------------------------
    # Final summary
    # ------------------------------------------------------------------

    print()

    print(
        f"Results directory: "
        f"{out_dir}"
    )

    print(
        f"Graphs discovered: "
        f"{len(graph_paths)}"
    )

    print(
        f"Graphs succeeded:  "
        f"{len(graph_paths) - len(failures)}"
    )

    print(
        f"Graphs failed:     "
        f"{len(failures)}"
    )

    print(
        f"Summary CSV:       "
        f"{summary_path}"
    )

    print(
        f"Summary rows:      "
        f"{len(all_rows)}"
    )

    if report_outputs is not None:

        print(
            f"Markdown report:   "
            f"{report_outputs['markdown']}"
        )

        if (
                report_outputs[
                    "html"
                ]
                is not None
        ):

            print(
                f"HTML report:       "
                f"{report_outputs['html']}"
            )

        else:

            print(
                "HTML report:       skipped"
            )

        if (
                report_outputs[
                    "pdf"
                ]
                is not None
        ):

            print(
                f"PDF report:        "
                f"{report_outputs['pdf']}"
            )

        else:

            print(
                "PDF report:        skipped"
            )

        print(
            f"Plots:             "
            f"{out_dir / 'plots'}"
        )

    elif report_error is not None:

        print(
            f"Report error:      "
            f"{report_error}"
        )

    # Graph experiment failures should still produce
    # a non-zero exit code.
    #
    # A skipped PDF export does NOT count as failure.

    if failures:
        raise SystemExit(1)

    if report_error is not None:
        raise SystemExit(1)


if __name__ == "__main__":
    main()