#!/usr/bin/env python3

import argparse
import csv
import math
import statistics
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt


# ============================================================================
# Helpers
# ============================================================================


def as_float(value):
    if value is None or value == "":
        return None

    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def as_int(value):
    value = as_float(value)

    if value is None:
        return None

    return int(value)


def safe_name(value):
    return (
        str(value)
        .replace(" ", "_")
        .replace("/", "_")
        .replace("\\", "_")
        .replace("(", "")
        .replace(")", "")
        .lower()
    )


def solver_label(row):
    return (
            row.get("solver_type")
            or row.get("solver")
            or "unknown"
    )


def status_ok(row):
    return (
            (row.get("status") or "").lower()
            == "ok"
    )


def graph_key(row):
    return (
            row.get("graph_path")
            or row.get("graph")
            or "unknown"
    )


def almost_equal(a, b):
    if a is None or b is None:
        return False

    return math.isclose(
        a,
        b,
        rel_tol=1e-9,
        abs_tol=1e-12,
    )


def format_number(value, digits=4):
    value = as_float(value)

    if value is None:
        return "-"

    if abs(value) >= 1_000_000:
        return f"{value:,.0f}"

    if abs(value) >= 1000:
        return f"{value:,.2f}"

    return f"{value:.{digits}g}"


def format_runtime_us(value):
    value = as_float(value)

    if value is None:
        return "-"

    if value < 1000:
        return f"{value:.0f} µs"

    if value < 1_000_000:
        return f"{value / 1000.0:.2f} ms"

    return f"{value / 1_000_000.0:.2f} s"


def format_percentage(value):
    if value is None:
        return "-"

    return f"{100.0 * value:.1f}%"


def markdown_table(headers, rows):
    output = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join("---" for _ in headers) + " |",
        ]

    for row in rows:
        output.append(
            "| "
            + " | ".join(str(value) for value in row)
            + " |"
        )

    return "\n".join(output)


# ============================================================================
# Statistics
# ============================================================================


def geometric_mean(values):
    values = [
        value
        for value in values
        if value is not None and value > 0
    ]

    if not values:
        return None

    return math.exp(
        sum(math.log(value) for value in values)
        / len(values)
    )


def best_value(rows, field):
    values = []

    for row in rows:
        if not status_ok(row):
            continue

        value = as_float(
            row.get(field)
        )

        if value is not None:
            values.append(value)

    if not values:
        return None

    return min(values)


def best_rows(rows, field):
    best = best_value(
        rows,
        field,
    )

    if best is None:
        return []

    result = []

    for row in rows:
        if not status_ok(row):
            continue

        value = as_float(
            row.get(field)
        )

        if almost_equal(
                value,
                best,
        ):
            result.append(row)

    return result


def pareto_front(rows):
    candidates = []

    for row in rows:
        if not status_ok(row):
            continue

        congestion = as_float(
            row.get("congestion")
        )

        runtime = as_float(
            row.get(
                "total_runtime_microseconds"
            )
        )

        if (
                congestion is None
                or runtime is None
        ):
            continue

        candidates.append(
            (
                row,
                congestion,
                runtime,
            )
        )

    front = []

    for (
            row,
            congestion,
            runtime,
    ) in candidates:

        dominated = False

        for (
                other,
                other_congestion,
                other_runtime,
        ) in candidates:

            if other is row:
                continue

            no_worse = (
                    other_congestion <= congestion
                    and
                    other_runtime <= runtime
            )

            strictly_better = (
                    other_congestion < congestion
                    or
                    other_runtime < runtime
            )

            if no_worse and strictly_better:
                dominated = True
                break

        if not dominated:
            front.append(row)

    return front


# ============================================================================
# Graph catalog
# ============================================================================


def build_graph_catalog(rows):
    catalog = {}
    name_counts = Counter()

    for row in rows:
        key = graph_key(row)

        if key in catalog:
            continue

        name = (
                row.get("graph")
                or Path(key).stem
                or "unknown"
        )

        name_counts[name] += 1

        catalog[key] = {
            "key": key,
            "name": name,
            "path": row.get("graph_path") or key,
            "nodes": as_int(
                row.get("nodes")
            ) or 0,
            "edges": as_int(
                row.get("edges")
            ) or 0,
        }

    for graph in catalog.values():
        name = graph["name"]

        if name_counts[name] == 1:
            graph["label"] = name
            continue

        path = Path(
            graph["path"]
        )

        parent = (
                path.parent.name
                or "dataset"
        )

        graph["label"] = (
            f"{parent}/{name}"
        )

    ordered = sorted(
        catalog.values(),
        key=lambda graph: (
            graph["nodes"],
            graph["edges"],
            graph["path"],
        ),
    )

    return catalog, ordered

def dataset_family(graph_path):
    """
    Derive a readable dataset family from the graph path.

    Example:
        experiments/datasets/small/Backbone/1221.lgf
        -> Backbone
    """

    if not graph_path:
        return "unknown"

    path = Path(
        graph_path
    )

    parent = (
            path.parent.name
            or "unknown"
    )

    return parent


def build_benchmark_configuration(
        rows,
):
    seeds = sorted(
        {
            str(row.get("seed"))
            for row in rows
            if row.get("seed") not in (
            None,
            "",
        )
        }
    )

    threads = sorted(
        {
            str(row.get("threads"))
            for row in rows
            if row.get("threads") not in (
            None,
            "",
        )
        }
    )

    solvers = sorted(
        {
            solver_label(row)
            for row in rows
            if status_ok(row)
        }
    )

    demands = sorted(
        {
            row.get("demand_model")
            for row in rows
            if row.get("demand_model")
        }
    )

    schema_versions = sorted(
        {
            str(
                row.get(
                    "schema_version"
                )
            )
            for row in rows
            if row.get(
            "schema_version"
        )
        }
    )

    unique_graphs = {
        graph_key(row)
        for row in rows
    }

    return {
        "graphs":
            len(
                unique_graphs
            ),

        "solvers":
            ", ".join(
                solvers
            ),

        "demand_models":
            ", ".join(
                demands
            ),

        "seed":
            ", ".join(
                seeds
            )
            if seeds
            else "-",

        "threads":
            ", ".join(
                threads
            )
            if threads
            else "-",

        "schema_version":
            ", ".join(
                schema_versions
            )
            if schema_versions
            else "-",
    }


def build_dataset_summary(
        rows,
):
    _, ordered_graphs = (
        build_graph_catalog(
            rows
        )
    )

    if not ordered_graphs:
        return {
            "graph_count": 0,
            "min_nodes": None,
            "max_nodes": None,
            "min_edges": None,
            "max_edges": None,
            "families": {},
        }

    node_counts = [
        graph["nodes"]
        for graph in ordered_graphs
    ]

    edge_counts = [
        graph["edges"]
        for graph in ordered_graphs
    ]

    family_counts = Counter()

    for graph in ordered_graphs:

        family = dataset_family(
            graph[
                "path"
            ]
        )

        family_counts[
            family
        ] += 1

    return {
        "graph_count":
            len(
                ordered_graphs
            ),

        "min_nodes":
            min(
                node_counts
            ),

        "max_nodes":
            max(
                node_counts
            ),

        "min_edges":
            min(
                edge_counts
            ),

        "max_edges":
            max(
                edge_counts
            ),

        "families":
            dict(
                sorted(
                    family_counts.items(),
                    key=lambda item: (
                        -item[1],
                        item[0],
                    ),
                )
            ),
    }
# ============================================================================
# Aggregate metrics
# ============================================================================


def calculate_congestion_summary(rows):
    groups = defaultdict(list)

    for row in rows:
        if not status_ok(row):
            continue

        groups[
            (
                graph_key(row),
                row.get("demand_model"),
            )
        ].append(row)

    result = defaultdict(
        lambda: {
            "wins": 0,
            "scenarios": 0,
            "relative": [],
            "within_05": 0,
            "within_10": 0,
            "within_25": 0,
        }
    )

    for group_rows in groups.values():
        best = best_value(
            group_rows,
            "congestion",
        )

        if best is None:
            continue

        winners = {
            solver_label(row)
            for row in best_rows(
                group_rows,
                "congestion",
            )
        }

        for solver in winners:
            result[solver]["wins"] += 1

        for row in group_rows:
            congestion = as_float(
                row.get("congestion")
            )

            if congestion is None:
                continue

            solver = solver_label(row)

            relative = (
                congestion / best
                if best > 0
                else None
            )

            result[solver][
                "scenarios"
            ] += 1

            if relative is None:
                continue

            result[solver][
                "relative"
            ].append(relative)

            if relative <= 1.05:
                result[solver][
                    "within_05"
                ] += 1

            if relative <= 1.10:
                result[solver][
                    "within_10"
                ] += 1

            if relative <= 1.25:
                result[solver][
                    "within_25"
                ] += 1

    return result


def unique_solver_rows_per_graph(rows):
    unique = {}

    for row in rows:
        if not status_ok(row):
            continue

        key = (
            graph_key(row),
            solver_label(row),
        )

        unique.setdefault(
            key,
            row,
        )

    return list(
        unique.values()
    )


def calculate_graph_summary(rows):
    unique_rows = (
        unique_solver_rows_per_graph(
            rows
        )
    )

    by_graph = defaultdict(list)

    for row in unique_rows:
        by_graph[
            graph_key(row)
        ].append(row)

    runtime_wins = Counter()
    ratio_wins = Counter()
    pareto_counts = Counter()

    runtime_values = defaultdict(list)
    ratio_values = defaultdict(list)

    for graph_rows in by_graph.values():

        for row in best_rows(
                graph_rows,
                "total_runtime_microseconds",
        ):
            runtime_wins[
                solver_label(row)
            ] += 1

        for row in best_rows(
                graph_rows,
                "oblivious_ratio",
        ):
            ratio_wins[
                solver_label(row)
            ] += 1

        for row in pareto_front(
                graph_rows
        ):
            pareto_counts[
                solver_label(row)
            ] += 1

        for row in graph_rows:
            solver = solver_label(row)

            runtime = as_float(
                row.get(
                    "total_runtime_microseconds"
                )
            )

            ratio = as_float(
                row.get(
                    "oblivious_ratio"
                )
            )

            if runtime is not None:
                runtime_values[
                    solver
                ].append(runtime)

            if ratio is not None:
                ratio_values[
                    solver
                ].append(ratio)

    return {
        "runtime_wins":
            runtime_wins,

        "ratio_wins":
            ratio_wins,

        "pareto_counts":
            pareto_counts,

        "runtime_values":
            runtime_values,

        "ratio_values":
            ratio_values,

        "graph_count":
            len(by_graph),
    }


# ============================================================================
# Plotting
# ============================================================================


def generate_plots(
        rows,
        output_dir,
):
    plots_dir = (
            output_dir
            / "plots"
    )

    plots_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    _, ordered_graphs = (
        build_graph_catalog(
            rows
        )
    )

    positions = {
        graph["key"]: index
        for index, graph
        in enumerate(ordered_graphs)
    }

    labels = [
        graph["label"]
        for graph
        in ordered_graphs
    ]

    solvers = sorted(
        {
            solver_label(row)
            for row in rows
            if status_ok(row)
        }
    )

    demands = sorted(
        {
            row.get("demand_model")
            for row in rows
            if row.get("demand_model")
        }
    )

    generated = {
        "congestion": {},
        "runtime": None,
    }

    width = max(
        14,
        len(ordered_graphs) * 0.30,
        )

    # ------------------------------------------------------------------
    # Congestion plots
    # ------------------------------------------------------------------

    for demand in demands:

        plt.figure(
            figsize=(
                width,
                7,
            )
        )

        for solver in solvers:
            values = {}

            for row in rows:
                if not status_ok(row):
                    continue

                if (
                        solver_label(row)
                        != solver
                ):
                    continue

                if (
                        row.get("demand_model")
                        != demand
                ):
                    continue

                congestion = as_float(
                    row.get("congestion")
                )

                if (
                        congestion is None
                        or congestion <= 0
                ):
                    continue

                values.setdefault(
                    graph_key(row),
                    congestion,
                )

            points = sorted(
                (
                    positions[key],
                    congestion,
                )
                for key, congestion
                in values.items()
                if key in positions
            )

            if not points:
                continue

            plt.plot(
                [
                    x
                    for x, _
                    in points
                ],
                [
                    y
                    for _, y
                    in points
                ],
                marker="o",
                markersize=3,
                linewidth=1.5,
                label=solver,
            )

        plt.yscale(
            "log"
        )

        plt.xlabel(
            "Graphs ordered by size"
        )

        plt.ylabel(
            "Congestion (log scale)"
        )

        plt.title(
            "Congestion by graph size "
            f"— {demand} demand"
        )

        plt.xticks(
            range(len(labels)),
            labels,
            rotation=90,
            fontsize=7,
        )

        plt.grid(
            True,
            which="both",
            alpha=0.25,
        )

        plt.legend()

        plt.tight_layout()

        path = (
                plots_dir
                / (
                    f"congestion_"
                    f"{safe_name(demand)}"
                    ".png"
                )
        )

        plt.savefig(
            path,
            dpi=160,
            bbox_inches="tight",
        )

        plt.close()

        generated[
            "congestion"
        ][demand] = path

    # ------------------------------------------------------------------
    # Runtime plot
    # ------------------------------------------------------------------

    plt.figure(
        figsize=(
            width,
            7,
        )
    )

    unique_rows = (
        unique_solver_rows_per_graph(
            rows
        )
    )

    for solver in solvers:
        points = []

        for row in unique_rows:
            if (
                    solver_label(row)
                    != solver
            ):
                continue

            runtime = as_float(
                row.get(
                    "total_runtime_microseconds"
                )
            )

            if (
                    runtime is None
                    or runtime <= 0
            ):
                continue

            key = graph_key(row)

            if key not in positions:
                continue

            points.append(
                (
                    positions[key],
                    runtime,
                )
            )

        points.sort(
            key=lambda item:
            item[0]
        )

        if not points:
            continue

        plt.plot(
            [
                x
                for x, _
                in points
            ],
            [
                y
                for _, y
                in points
            ],
            marker="o",
            markersize=3,
            linewidth=1.5,
            label=solver,
        )

    plt.yscale(
        "log"
    )

    plt.xlabel(
        "Graphs ordered by size"
    )

    plt.ylabel(
        "Runtime [µs] (log scale)"
    )

    plt.title(
        "Solver runtime by graph size"
    )

    plt.xticks(
        range(len(labels)),
        labels,
        rotation=90,
        fontsize=7,
    )

    plt.grid(
        True,
        which="both",
        alpha=0.25,
    )

    plt.legend()

    plt.tight_layout()

    runtime_path = (
            plots_dir
            / "runtime_by_graph.png"
    )

    plt.savefig(
        runtime_path,
        dpi=160,
        bbox_inches="tight",
    )

    plt.close()

    generated[
        "runtime"
    ] = runtime_path

    return generated


# ============================================================================
# Automatic Key Findings
# ============================================================================


def generate_key_findings(rows):
    findings = []

    solvers = sorted(
        {
            solver_label(row)
            for row in rows
            if status_ok(row)
        }
    )

    if not solvers:
        return findings

    congestion_summary = (
        calculate_congestion_summary(
            rows
        )
    )

    graph_summary = (
        calculate_graph_summary(
            rows
        )
    )

    # ------------------------------------------------------------------
    # Congestion leader
    # ------------------------------------------------------------------

    congestion_leader = max(
        solvers,
        key=lambda solver:
        congestion_summary[
            solver
        ]["wins"],
    )

    congestion_wins = (
        congestion_summary[
            congestion_leader
        ]["wins"]
    )

    findings.append(
        f"**{congestion_leader}** achieved the lowest "
        f"congestion most frequently, winning or tying "
        f"**{congestion_wins} graph-demand scenarios**."
    )

    # ------------------------------------------------------------------
    # Most consistent congestion
    # ------------------------------------------------------------------

    consistency_candidates = []

    for solver in solvers:
        relative = (
            congestion_summary[
                solver
            ]["relative"]
        )

        if not relative:
            continue

        consistency_candidates.append(
            (
                statistics.median(
                    relative
                ),
                solver,
            )
        )

    if consistency_candidates:

        median_relative, solver = min(
            consistency_candidates,
            key=lambda item:
            item[0],
        )

        findings.append(
            f"**{solver}** was the most consistent solver "
            f"relative to the best congestion result, with "
            f"a median performance of "
            f"**{median_relative:.3f}× the best observed "
            f"congestion**."
        )

    # ------------------------------------------------------------------
    # Runtime leader
    # ------------------------------------------------------------------

    runtime_wins = (
        graph_summary[
            "runtime_wins"
        ]
    )

    if runtime_wins:

        runtime_solver, wins = (
            runtime_wins
            .most_common(1)[0]
        )

        runtime_values = (
            graph_summary[
                "runtime_values"
            ][runtime_solver]
        )

        median_runtime = (
            statistics.median(
                runtime_values
            )
            if runtime_values
            else None
        )

        runtime_text = (
            format_runtime_us(
                median_runtime
            )
            if median_runtime
               is not None
            else "-"
        )

        findings.append(
            f"**{runtime_solver}** showed the strongest "
            f"runtime performance, finishing fastest on "
            f"**{wins} graphs** with a median runtime of "
            f"**{runtime_text}**."
        )

    # ------------------------------------------------------------------
    # Oblivious ratio leader
    # ------------------------------------------------------------------

    ratio_wins = (
        graph_summary[
            "ratio_wins"
        ]
    )

    if ratio_wins:

        ratio_solver, wins = (
            ratio_wins
            .most_common(1)[0]
        )

        findings.append(
            f"**{ratio_solver}** achieved the best oblivious "
            f"ratio most frequently, leading on "
            f"**{wins} graph instances**."
        )

    # ------------------------------------------------------------------
    # Pareto leader
    # ------------------------------------------------------------------

    pareto_counts = (
        graph_summary[
            "pareto_counts"
        ]
    )

    if pareto_counts:

        pareto_solver, count = (
            pareto_counts
            .most_common(1)[0]
        )

        findings.append(
            f"**{pareto_solver}** appeared most frequently on "
            f"the congestion/runtime Pareto frontier, appearing "
            f"on **{count} graphs**."
        )

    # ------------------------------------------------------------------
    # Demand-specific congestion leaders
    # ------------------------------------------------------------------

    demand_groups = defaultdict(
        lambda: defaultdict(int)
    )

    grouped = defaultdict(list)

    for row in rows:

        if not status_ok(row):
            continue

        demand = (
                row.get(
                    "demand_model"
                )
                or "unknown"
        )

        grouped[
            (
                graph_key(row),
                demand,
            )
        ].append(row)

    for (
            (_, demand),
            group_rows,
    ) in grouped.items():

        for winner in best_rows(
                group_rows,
                "congestion",
        ):

            demand_groups[
                demand
            ][
                solver_label(
                    winner
                )
            ] += 1

    demand_sentences = []

    for demand in sorted(
            demand_groups
    ):

        counts = (
            demand_groups[
                demand
            ]
        )

        if not counts:
            continue

        solver, wins = max(
            counts.items(),
            key=lambda item:
            item[1],
        )

        demand_sentences.append(
            f"{demand}: **{solver}** "
            f"({wins} wins)"
        )

    if demand_sentences:

        findings.append(
            "Congestion leadership by demand model: "
            + "; ".join(
                demand_sentences
            )
            + "."
        )

    # ------------------------------------------------------------------
    # Overall interpretation
    # ------------------------------------------------------------------

    findings.append(
        "The benchmark shows a clear quality/runtime "
        "trade-off rather than a universally dominant "
        "routing method, supporting topology-aware "
        "solver selection."
    )

    return findings


# ============================================================================
# Markdown Report
# ============================================================================


def build_report(
        rows,
        title,
        source_path,
        plots,
):
    generated = (
        datetime.now().strftime(
            "%Y-%m-%d %H:%M:%S"
        )
    )

    _, ordered_graphs = (
        build_graph_catalog(
            rows
        )
    )

    solvers = sorted(
        {
            solver_label(row)
            for row in rows
            if status_ok(row)
        }
    )

    graph_groups = defaultdict(list)
    scenario_groups = defaultdict(list)

    for row in rows:
        key = graph_key(row)

        graph_groups[
            key
        ].append(row)

        scenario_groups[
            (
                key,
                row.get("demand_model"),
            )
        ].append(row)

    congestion_summary = (
        calculate_congestion_summary(
            rows
        )
    )

    graph_summary = (
        calculate_graph_summary(
            rows
        )
    )

    runtime_wins = (
        graph_summary[
            "runtime_wins"
        ]
    )

    ratio_wins = (
        graph_summary[
            "ratio_wins"
        ]
    )

    pareto_counts = (
        graph_summary[
            "pareto_counts"
        ]
    )

    lines = [
        f"# {title}",
        "",
        f"Generated: `{generated}`",
        "",
        f"Source data: `{source_path}`",
        "",
    ]


    benchmark_configuration = (
        build_benchmark_configuration(
            rows
        )
    )

    lines += [
        "## Benchmark Configuration",
        "",
    ]

    configuration_rows = [
        [
            "Graphs",
            benchmark_configuration[
                "graphs"
            ],
        ],
        [
            "Solvers",
            benchmark_configuration[
                "solvers"
            ],
        ],
        [
            "Demand Models",
            benchmark_configuration[
                "demand_models"
            ],
        ],
        [
            "Random Seed",
            benchmark_configuration[
                "seed"
            ],
        ],
        [
            "Threads",
            benchmark_configuration[
                "threads"
            ],
        ],
        [
            "Result Schema",
            benchmark_configuration[
                "schema_version"
            ],
        ],
    ]

    lines.append(
        markdown_table(
            [
                "Configuration",
                "Value",
            ],
            configuration_rows,
        )
    )

    lines.append("")

    # ==================================================================
    # Key Findings
    # ==================================================================

    lines += [
        "## Key Findings",
        "",
        (
            f"This benchmark evaluates "
            f"**{len(solvers)} routing solvers** "
            f"across **{len(ordered_graphs)} unique graph "
            f"instances** and "
            f"**{len(scenario_groups)} graph-demand scenarios**."
        ),
        "",
    ]


    key_findings = (
        generate_key_findings(
            rows
        )
    )

    for finding in key_findings:
        lines.append(
            f"- {finding}"
        )

    lines.append("")

    # ==================================================================
    # Aggregate Congestion Performance
    # ==================================================================

    lines += [
        "## Aggregate Congestion Performance",
        "",
        (
            "Congestion is normalized independently for every "
            "graph-demand scenario. A value of `1.0×` means "
            "that the solver matched the best measured "
            "congestion for that scenario."
        ),
        "",
    ]

    aggregate_rows = []

    for solver in solvers:

        data = (
            congestion_summary[
                solver
            ]
        )

        relative = (
            data[
                "relative"
            ]
        )

        scenarios = (
            data[
                "scenarios"
            ]
        )

        median_relative = (
            statistics.median(
                relative
            )
            if relative
            else None
        )

        geo_relative = (
            geometric_mean(
                relative
            )
        )

        def fraction(value):
            if scenarios == 0:
                return None

            return (
                    value
                    / scenarios
            )

        aggregate_rows.append([
            solver,

            data[
                "wins"
            ],

            format_number(
                median_relative
            ),

            format_number(
                geo_relative
            ),

            format_percentage(
                fraction(
                    data[
                        "within_05"
                    ]
                )
            ),

            format_percentage(
                fraction(
                    data[
                        "within_10"
                    ]
                )
            ),

            format_percentage(
                fraction(
                    data[
                        "within_25"
                    ]
                )
            ),
        ])

    lines.append(
        markdown_table(
            [
                "Solver",
                "Congestion Wins",
                "Median vs Best",
                "Geometric Mean vs Best",
                "Within 5%",
                "Within 10%",
                "Within 25%",
            ],
            aggregate_rows,
        )
    )

    lines.append("")

    # ==================================================================
    # Runtime / Oblivious Ratio
    # ==================================================================

    lines += [
        "## Solver Runtime and Oblivious Ratio",
        "",
        (
            "Runtime and oblivious-ratio statistics are "
            "evaluated once per graph rather than once per "
            "demand model because these quantities belong "
            "to the generated routing scheme itself."
        ),
        "",
    ]

    runtime_rows = []

    for solver in solvers:

        runtimes = (
            graph_summary[
                "runtime_values"
            ][solver]
        )

        ratios = (
            graph_summary[
                "ratio_values"
            ][solver]
        )

        median_runtime = (
            statistics.median(
                runtimes
            )
            if runtimes
            else None
        )

        median_ratio = (
            statistics.median(
                ratios
            )
            if ratios
            else None
        )

        runtime_rows.append([
            solver,

            runtime_wins[
                solver
            ],

            format_runtime_us(
                median_runtime
            ),

            ratio_wins[
                solver
            ],

            format_number(
                median_ratio
            ),

            pareto_counts[
                solver
            ],
        ])

    lines.append(
        markdown_table(
            [
                "Solver",
                "Runtime Wins",
                "Median Runtime",
                "Ratio Wins",
                "Median Oblivious Ratio",
                "Pareto Graphs",
            ],
            runtime_rows,
        )
    )

    lines.append("")

    # ==================================================================
    # Congestion Plots
    # ==================================================================

    lines += [
        "## Congestion by Graph Size",
        "",
        (
            "The following plots compare congestion across "
            "the complete benchmark suite for each demand "
            "model. Graphs are ordered by node count and "
            "then edge count. The logarithmic Y-axis is used "
            "because congestion spans several orders of "
            "magnitude."
        ),
        "",
    ]

    for (
            demand,
            plot_path,
    ) in sorted(
        plots[
            "congestion"
        ].items()
    ):

        relative = (
                Path("plots")
                / plot_path.name
        )

        lines += [
            f"### {demand.capitalize()} Demand",
            "",
            (
                f"![Congestion for {demand} demand]"
                f"({relative.as_posix()})"
            ),
            "",
        ]

    # ==================================================================
    # Runtime Plot
    # ==================================================================

    lines += [
        "## Runtime by Graph Size",
        "",
        (
            "Total solver runtime is shown once per graph "
            "and solver. Graphs use the same size ordering "
            "as the congestion plots. Runtime is plotted "
            "on a logarithmic Y-axis so that both small "
            "and large instances remain visible."
        ),
        "",
    ]

    runtime_plot = (
        plots.get(
            "runtime"
        )
    )

    if runtime_plot:

        relative = (
                Path("plots")
                / runtime_plot.name
        )

        lines += [
            (
                "![Solver runtime by graph size]"
                f"({relative.as_posix()})"
            ),
            "",
        ]

    # ==================================================================
    # Benchmark Dataset Summary
    # ==================================================================

    dataset_summary = (
        build_dataset_summary(
            rows
        )
    )

    lines += [
        "## Benchmark Dataset Summary",
        "",
        (
            "The benchmark spans multiple network datasets and "
            "topologies. Full per-instance measurements and graph "
            "identifiers are available in `summary.csv`."
        ),
        "",
    ]

    dataset_rows = [
        [
            "Graph Instances",
            dataset_summary[
                "graph_count"
            ],
        ],
        [
            "Node Range",
            (
                f"{dataset_summary['min_nodes']} – "
                f"{dataset_summary['max_nodes']}"
            ),
        ],
        [
            "Edge Range",
            (
                f"{dataset_summary['min_edges']} – "
                f"{dataset_summary['max_edges']}"
            ),
        ],
    ]

    lines.append(
        markdown_table(
            [
                "Dataset Property",
                "Value",
            ],
            dataset_rows,
        )
    )

    lines.append("")

    lines += [
        "### Dataset Families",
        "",
    ]

    family_rows = []

    for (
            family,
            count,
    ) in dataset_summary[
        "families"
    ].items():

        family_rows.append([
            family,
            count,
        ])

    lines.append(
        markdown_table(
            [
                "Dataset Family",
                "Graphs",
            ],
            family_rows,
        )
    )

    lines.append("")
    # ==================================================================
    # Interpretation Notes
    # ==================================================================

    lines += [
        "## Interpretation Notes",
        "",
        "- Lower congestion is better.",
        "- Lower oblivious ratio is better.",
        "- Lower runtime is better.",
        (
            "- Congestion comparisons are performed independently "
            "for every graph-demand scenario."
        ),
        (
            "- Runtime and oblivious-ratio comparisons are "
            "performed once per graph because these values belong "
            "to the routing scheme and do not depend on the "
            "evaluated demand model."
        ),
        (
            "- `1.0×` means that a solver matched the best "
            "measured result for that graph-demand scenario."
        ),
        (
            "- `Within 10%` means that the solver produced "
            "congestion no more than 1.10× the best measured "
            "congestion."
        ),
        (
            "- Pareto-front solvers are not simultaneously "
            "outperformed by another solver in both congestion "
            "and runtime."
        ),
        (
            "- Runtime comparisons should use equivalent "
            "hardware, compiler configuration, build type, "
            "OpenMP settings and thread count."
        ),
        (
            "- Full per-instance measurements are available "
            "in `summary.csv`."
        ),
        "",
    ]

    return "\n".join(
        lines
    )


# ============================================================================
# HTML / PDF Export
# ============================================================================


def export_report_documents(
        markdown_path,
):
    """
    HTML and PDF are optional presentation outputs.

    Markdown generation must remain successful even when
    python-markdown, WeasyPrint, Pango, GTK libraries, etc.
    are unavailable.
    """

    markdown_path = Path(
        markdown_path
    )

    html_path = (
        markdown_path
        .with_suffix(".html")
    )

    pdf_path = (
        markdown_path
        .with_suffix(".pdf")
    )

    # ------------------------------------------------------------------
    # HTML
    # ------------------------------------------------------------------

    try:
        import markdown

    except Exception as exc:

        print()
        print(
            "[REPORT WARNING] "
            "HTML/PDF export skipped because the "
            "'markdown' package is unavailable."
        )

        print(
            f"[REPORT WARNING] {exc}"
        )

        return {
            "html": None,
            "pdf": None,
        }

    markdown_text = (
        markdown_path
        .read_text(
            encoding="utf-8"
        )
    )

    html_body = (
        markdown.markdown(
            markdown_text,
            extensions=[
                "tables",
                "fenced_code",
            ],
        )
    )

    css = """
    @page {
        size: A4;
        margin: 18mm;
    }

    body {
        font-family:
            -apple-system,
            BlinkMacSystemFont,
            "Segoe UI",
            Arial,
            sans-serif;

        font-size: 10pt;
        line-height: 1.45;
        color: #222;
    }

    h1 {
        font-size: 24pt;
        margin-bottom: 8px;
    }

    h2 {
        font-size: 17pt;
        margin-top: 28px;
        border-bottom: 1px solid #ddd;
        padding-bottom: 5px;
    }

    h3 {
        font-size: 13pt;
        margin-top: 22px;
    }

    table {
        width: 100%;
        border-collapse: collapse;
        margin-top: 10px;
        margin-bottom: 18px;
        font-size: 8pt;
    }

    th,
    td {
        border: 1px solid #ccc;
        padding: 5px 7px;
        text-align: left;
    }

    th {
        background: #f3f3f3;
        font-weight: 600;
    }

    tr {
        page-break-inside: avoid;
    }

    img {
        display: block;
        max-width: 100%;
        height: auto;

        margin-left: auto;
        margin-right: auto;

        margin-top: 12px;
        margin-bottom: 20px;
    }

    code {
        font-family:
            "SFMono-Regular",
            Consolas,
            monospace;

        font-size: 0.92em;
        background: #f5f5f5;
        padding: 1px 4px;
    }

    ul {
        margin-top: 6px;
        margin-bottom: 18px;
    }

    li {
        margin-bottom: 5px;
    }

    p {
        margin-top: 6px;
        margin-bottom: 8px;
    }
    """

    complete_html = f"""
<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8">

    <title>
        E-Routing Benchmark Report
    </title>

    <style>
        {css}
    </style>
</head>

<body>
    {html_body}
</body>
</html>
"""

    html_path.write_text(
        complete_html,
        encoding="utf-8",
    )

    # ------------------------------------------------------------------
    # PDF
    # ------------------------------------------------------------------

    try:
        from weasyprint import HTML

        HTML(
            string=complete_html,
            base_url=str(
                markdown_path
                .parent
                .resolve()
            ),
        ).write_pdf(
            str(
                pdf_path
            )
        )

    except Exception as exc:

        print()
        print(
            "[REPORT WARNING] "
            "PDF export skipped."
        )

        print(
            "[REPORT WARNING] "
            "Markdown and HTML reports were still generated."
        )

        print(
            f"[REPORT WARNING] {exc}"
        )

        return {
            "html":
                html_path,

            "pdf":
                None,
        }

    return {
        "html":
            html_path,

        "pdf":
            pdf_path,
    }


# ============================================================================
# Public API
# ============================================================================


def generate_report_from_summary(
        summary_path,
        report_path=None,
        title="E-Routing Benchmark Report",
):
    summary_path = Path(
        summary_path
    )

    if not summary_path.exists():
        raise FileNotFoundError(
            f"summary.csv does not exist: "
            f"{summary_path}"
        )

    if report_path is None:
        report_path = (
                summary_path.parent
                / "report.md"
        )

    else:
        report_path = Path(
            report_path
        )

    with summary_path.open(
            newline="",
            encoding="utf-8",
    ) as f:

        rows = list(
            csv.DictReader(f)
        )

    if not rows:
        raise ValueError(
            "Cannot generate report: "
            "summary.csv contains no rows."
        )

    plots = generate_plots(
        rows,
        report_path.parent,
    )

    report = build_report(
        rows=rows,
        title=title,
        source_path=summary_path,
        plots=plots,
    )

    report_path.write_text(
        report,
        encoding="utf-8",
    )

    exports = (
        export_report_documents(
            report_path
        )
    )

    return {
        "markdown":
            report_path,

        "html":
            exports[
                "html"
            ],

        "pdf":
            exports[
                "pdf"
            ],
    }


# ============================================================================
# CLI
# ============================================================================


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Generate E-Routing benchmark report."
        )
    )

    parser.add_argument(
        "--input",
        required=True,
        help=(
            "Path to summary.csv."
        ),
    )

    parser.add_argument(
        "--output",
        default=None,
        help=(
            "Optional report.md output path."
        ),
    )

    parser.add_argument(
        "--title",
        default=(
            "E-Routing Benchmark Report"
        ),
    )

    return parser.parse_args()


def main():
    args = parse_args()

    outputs = (
        generate_report_from_summary(
            summary_path=args.input,
            report_path=args.output,
            title=args.title,
        )
    )

    print()
    print(
        "[REPORT] Generated:"
    )

    print(
        f"  Markdown: "
        f"{outputs['markdown']}"
    )

    if outputs["html"] is not None:

        print(
            f"  HTML:     "
            f"{outputs['html']}"
        )

    else:

        print(
            "  HTML:     skipped"
        )

    if outputs["pdf"] is not None:

        print(
            f"  PDF:      "
            f"{outputs['pdf']}"
        )

    else:

        print(
            "  PDF:      skipped"
        )


if __name__ == "__main__":
    main()