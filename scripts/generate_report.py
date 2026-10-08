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

def as_positive_float(value):
    value = as_float(value)

    if value is None or value <= 0:
        return None

    return value

def as_int(value):
    value = as_float(value)
    return None if value is None else int(value)


def safe_name(value):
    return (
        str(value)
        .replace(" ", "_")
        .replace("/", "_")
        .replace("\\", "_")
        .replace("(", "")
        .replace(")", "")
        .replace(":", "_")
        .lower()
    )


def solver_label(row):
    return row.get("solver_type") or row.get("solver") or "unknown"


def status_ok(row):
    return (row.get("status") or "").lower() == "ok"


def graph_key(row):
    return row.get("graph_path") or row.get("graph") or "unknown"


def almost_equal(a, b):
    if a is None or b is None:
        return False
    return math.isclose(a, b, rel_tol=1e-9, abs_tol=1e-12)

def format_positive_number(
        value,
        digits=4,
):
    value = as_positive_float(
        value
    )

    if value is None:
        return "N/A"

    return format_number(
        value,
        digits,
    )

def has_valid_failure_analysis(row):
    tested_links = as_int(
        row.get(
            "failure_tested_links"
        )
    )

    return (
            status_ok(row)
            and
            tested_links is not None
            and
            tested_links > 0
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
    if value is None or value < 0:
        return "N/A"
    if value < 1000:
        return f"{value:.0f} us"
    if value < 1_000_000:
        return f"{value / 1000.0:.2f} ms"
    return f"{value / 1_000_000.0:.2f} s"


def format_percentage(value):
    value = as_float(value)
    if value is None or value < 0:
        return "N/A"
    return f"{100.0 * value:.1f}%"


def format_edge(source, target):
    source = as_int(source)
    target = as_int(target)

    if source is None or target is None or source < 0 or target < 0:
        return "N/A"

    return f"{source} -> {target}"


def median_or_none(values):
    values = [v for v in values if v is not None]
    return statistics.median(values) if values else None

def has_valid_failure_analysis(row):
    tested_links = as_int(
        row.get(
            "failure_tested_links"
        )
    )

    return (
            status_ok(row)
            and
            tested_links is not None
            and
            tested_links > 0
    )

def has_valid_failure_recovery(row):
    tested_links = as_int(
        row.get(
            "recovery_tested_links"
        )
    )

    return (
            status_ok(row)
            and
            tested_links is not None
            and
            tested_links > 0
    )


def has_recovery_data(rows):
    return any(
        has_valid_failure_recovery(row)
        for row in rows
    )


def geometric_mean(values):
    values = [v for v in values if v is not None and v > 0]
    if not values:
        return None
    return math.exp(sum(math.log(v) for v in values) / len(values))


def markdown_table(headers, rows):
    out = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join("---" for _ in headers) + " |",
        ]
    for row in rows:
        out.append("| " + " | ".join(str(v) for v in row) + " |")
    return "\n".join(out)


def best_value(rows, field):
    vals = [
        as_float(row.get(field))
        for row in rows
        if status_ok(row) and as_float(row.get(field)) is not None
    ]
    return min(vals) if vals else None


def best_rows(rows, field):
    best = best_value(rows, field)
    if best is None:
        return []
    return [
        row for row in rows
        if status_ok(row)
           and almost_equal(as_float(row.get(field)), best)
    ]


# ============================================================================
# Graph catalog / benchmark configuration
# ============================================================================

def build_graph_catalog(rows):
    catalog = {}
    name_counts = Counter()

    for row in rows:
        key = graph_key(row)
        if key in catalog:
            continue

        name = row.get("graph") or Path(key).stem or "unknown"
        name_counts[name] += 1
        catalog[key] = {
            "key": key,
            "name": name,
            "path": row.get("graph_path") or key,
            "nodes": as_int(row.get("nodes")) or 0,
            "edges": as_int(row.get("edges")) or 0,
        }

    for graph in catalog.values():
        name = graph["name"]
        if name_counts[name] == 1:
            graph["label"] = name
        else:
            graph["label"] = f"{Path(graph['path']).parent.name}/{name}"

    ordered = sorted(
        catalog.values(),
        key=lambda g: (g["nodes"], g["edges"], g["path"]),
    )
    return catalog, ordered


def dataset_family(graph_path):
    if not graph_path:
        return "unknown"
    return Path(graph_path).parent.name or "unknown"


def build_benchmark_configuration(rows):
    return {
        "graphs": len({graph_key(row) for row in rows}),
        "solvers": ", ".join(sorted({
            solver_label(row) for row in rows if status_ok(row)
        })),
        "demand_models": ", ".join(sorted({
            row.get("demand_model")
            for row in rows if row.get("demand_model")
        })),
        "seed": ", ".join(sorted({
            str(row.get("seed"))
            for row in rows if row.get("seed") not in (None, "")
        })) or "-",
        "threads": ", ".join(sorted({
            str(row.get("threads"))
            for row in rows if row.get("threads") not in (None, "")
        })) or "-",
        "schema_version": ", ".join(sorted({
            str(row.get("schema_version"))
            for row in rows if row.get("schema_version")
        })) or "-",
        "failure_recovery": ", ".join(sorted({
            str(row.get("failure_recovery"))
            for row in rows
            if row.get("failure_recovery") not in (None, "")
        })) or "-",
    }


def build_dataset_summary(rows):
    _, graphs = build_graph_catalog(rows)
    if not graphs:
        return {
            "graph_count": 0,
            "min_nodes": None,
            "max_nodes": None,
            "min_edges": None,
            "max_edges": None,
            "families": {},
        }

    family_counts = Counter(dataset_family(g["path"]) for g in graphs)
    return {
        "graph_count": len(graphs),
        "min_nodes": min(g["nodes"] for g in graphs),
        "max_nodes": max(g["nodes"] for g in graphs),
        "min_edges": min(g["edges"] for g in graphs),
        "max_edges": max(g["edges"] for g in graphs),
        "families": dict(sorted(
            family_counts.items(),
            key=lambda item: (-item[1], item[0]),
        )),
    }


# ============================================================================
# Benchmark statistics
# ============================================================================

def calculate_congestion_summary(rows):
    groups = defaultdict(list)
    for row in rows:
        if status_ok(row):
            groups[(graph_key(row), row.get("demand_model"))].append(row)

    result = defaultdict(lambda: {
        "wins": 0,
        "scenarios": 0,
        "relative": [],
        "within_05": 0,
        "within_10": 0,
        "within_25": 0,
    })

    for group_rows in groups.values():
        best = best_value(group_rows, "congestion")
        if best is None:
            continue

        for solver in {
            solver_label(row)
            for row in best_rows(group_rows, "congestion")
        }:
            result[solver]["wins"] += 1

        for row in group_rows:
            congestion = as_float(row.get("congestion"))
            if congestion is None:
                continue

            solver = solver_label(row)
            result[solver]["scenarios"] += 1

            relative = congestion / best if best > 0 else None
            if relative is None:
                continue

            result[solver]["relative"].append(relative)
            if relative <= 1.05:
                result[solver]["within_05"] += 1
            if relative <= 1.10:
                result[solver]["within_10"] += 1
            if relative <= 1.25:
                result[solver]["within_25"] += 1

    return result


def unique_solver_rows_per_graph(rows):
    unique = {}
    for row in rows:
        if status_ok(row):
            unique.setdefault((graph_key(row), solver_label(row)), row)
    return list(unique.values())


def calculate_graph_summary(rows):
    unique_rows = unique_solver_rows_per_graph(rows)
    by_graph = defaultdict(list)
    for row in unique_rows:
        by_graph[graph_key(row)].append(row)

    runtime_wins = Counter()
    ratio_wins = Counter()
    runtime_values = defaultdict(list)
    ratio_values = defaultdict(list)

    for graph_rows in by_graph.values():
        for row in best_rows(graph_rows, "total_runtime_microseconds"):
            runtime_wins[solver_label(row)] += 1
            valid_ratio_rows = [row
                                for row in graph_rows
                                if (
                                        as_positive_float(
                                            row.get(
                                                "oblivious_ratio"
                                            )
                                        )
                                        is not None
                                )
                                ]

            for row in best_rows(
                    valid_ratio_rows,
                    "oblivious_ratio",
            ):
                ratio_wins[
                    solver_label(row)
                ] += 1

        for row in graph_rows:
            solver = solver_label(row)
            runtime = as_float(row.get("total_runtime_microseconds"))
            ratio = as_positive_float(
                row.get(
                    "oblivious_ratio"
                )
            )
            if runtime is not None:
                runtime_values[solver].append(runtime)
            if ratio is not None:
                ratio_values[solver].append(ratio)

    return {
        "runtime_wins": runtime_wins,
        "ratio_wins": ratio_wins,
        "runtime_values": runtime_values,
        "ratio_values": ratio_values,
    }


# ============================================================================
# Failure / resilience statistics
# ============================================================================

def has_failure_data(rows):
    return any(
        has_valid_failure_analysis(
            row
        )
        for row in rows
    )


def calculate_failure_summary(rows):
    scenario_groups = defaultdict(list)

    for row in rows:

        if not has_valid_failure_analysis(row):
            continue

        scenario_groups[
            (
                graph_key(row),
                row.get(
                    "demand_model"
                ),
            )
        ].append(row)

    summary = defaultdict(lambda: {
        "wins": 0,
        "scenarios": 0,
        "max_loss": [],
        "avg_loss": [],
        "median_loss": [],
        "max_affected": [],
        "critical_10": [],
        "critical_25": [],
        "critical_50": [],
        "failure_runtime": [],
        "evaluation_runtime": [],
    })

    for scenario_rows in scenario_groups.values():
        for row in best_rows(
                scenario_rows,
                "failure_maximum_lost_traffic_fraction",
        ):
            summary[solver_label(row)]["wins"] += 1

        for row in scenario_rows:
            solver = solver_label(row)
            data = summary[solver]
            data["scenarios"] += 1

            mappings = {
                "max_loss": "failure_maximum_lost_traffic_fraction",
                "avg_loss": "failure_average_lost_traffic_fraction",
                "median_loss": "failure_median_lost_traffic_fraction",
                "max_affected": "failure_maximum_affected_demand_fraction",
                "critical_10": "failure_critical_links_10_percent",
                "critical_25": "failure_critical_links_25_percent",
                "critical_50": "failure_critical_links_50_percent",
                "failure_runtime": "failure_analysis_runtime_microseconds",
                "evaluation_runtime": "evaluation_runtime_microseconds",
            }

            for target, field in mappings.items():
                value = as_float(row.get(field))
                if value is not None:
                    data[target].append(value)

    return summary


def top_critical_scenarios(rows, limit=15):
    candidates = []
    for row in rows:

        if not has_valid_failure_analysis(row):
            continue
        loss = as_float(row.get("failure_maximum_lost_traffic_fraction"))
        if loss is None:
            continue
        candidates.append((loss, row))

    candidates.sort(key=lambda item: item[0], reverse=True)
    return [row for _, row in candidates[:limit]]



def calculate_recovery_summary(rows):
    scenario_groups = defaultdict(list)

    for row in rows:
        if not has_valid_failure_recovery(row):
            continue

        scenario_groups[
            (
                graph_key(row),
                row.get("demand_model"),
            )
        ].append(row)

    summary = defaultdict(lambda: {
        "quality_wins": 0,
        "speed_wins": 0,
        "scenarios": 0,
        "relative_worst_post": [],
        "max_post": [],
        "avg_post": [],
        "max_factor": [],
        "avg_factor": [],
        "max_unroutable": [],
        "avg_unroutable": [],
        "disconnect_rate": [],
        "recoverability_rate": [],
        "avg_runtime": [],
        "max_runtime": [],
    })

    for scenario_rows in scenario_groups.values():
        quality_rows = [
            row
            for row in scenario_rows
            if (
                    as_positive_float(
                        row.get(
                            "recovery_maximum_post_failure_congestion"
                        )
                    )
                    is not None
            )
        ]

        if quality_rows:
            best_quality = best_value(
                quality_rows,
                "recovery_maximum_post_failure_congestion",
            )

            for row in best_rows(
                    quality_rows,
                    "recovery_maximum_post_failure_congestion",
            ):
                summary[
                    solver_label(row)
                ]["quality_wins"] += 1
        else:
            best_quality = None

        speed_rows = [
            row
            for row in scenario_rows
            if (
                    as_positive_float(
                        row.get(
                            "recovery_average_recomputation_runtime_microseconds"
                        )
                    )
                    is not None
            )
        ]

        for row in best_rows(
                speed_rows,
                "recovery_average_recomputation_runtime_microseconds",
        ):
            summary[
                solver_label(row)
            ]["speed_wins"] += 1

        for row in scenario_rows:
            solver = solver_label(row)
            data = summary[solver]
            data["scenarios"] += 1

            max_post = as_positive_float(
                row.get(
                    "recovery_maximum_post_failure_congestion"
                )
            )

            if (
                    max_post is not None
                    and best_quality is not None
                    and best_quality > 0
            ):
                data["relative_worst_post"].append(
                    max_post / best_quality
                )

            mappings = {
                "max_post":
                    "recovery_maximum_post_failure_congestion",
                "avg_post":
                    "recovery_average_post_failure_congestion",
                "max_factor":
                    "recovery_maximum_congestion_increase_factor",
                "avg_factor":
                    "recovery_average_congestion_increase_factor",
                "max_unroutable":
                    "recovery_maximum_unroutable_demand_fraction",
                "avg_unroutable":
                    "recovery_average_unroutable_demand_fraction",
                "avg_runtime":
                    "recovery_average_recomputation_runtime_microseconds",
                "max_runtime":
                    "recovery_maximum_recomputation_runtime_microseconds",
            }

            for target, field in mappings.items():
                value = as_float(
                    row.get(field)
                )
                if value is not None and value >= 0:
                    data[target].append(value)

            tested = as_float(
                row.get(
                    "recovery_tested_links"
                )
            )
            disconnected = as_float(
                row.get(
                    "recovery_disconnected_failures"
                )
            )
            successful = as_float(
                row.get(
                    "recovery_successful_recomputations"
                )
            )

            if (
                    tested is not None
                    and tested > 0
                    and disconnected is not None
            ):
                data["disconnect_rate"].append(
                    disconnected / tested
                )

                survivable = (
                        tested
                        -
                        disconnected
                )

                if (
                        survivable > 0
                        and successful is not None
                ):
                    data["recoverability_rate"].append(
                        successful
                        /
                        survivable
                    )

    return summary


def top_recovery_scenarios(rows, limit=15):
    candidates = []

    for row in rows:
        if not has_valid_failure_recovery(row):
            continue

        degradation = as_float(
            row.get(
                "recovery_maximum_congestion_increase_factor"
            )
        )

        if (
                degradation is None
                or degradation < 0
        ):
            continue

        candidates.append(
            (
                degradation,
                row,
            )
        )

    candidates.sort(
        key=lambda item: item[0],
        reverse=True,
    )

    return [
        row
        for _, row in candidates[:limit]
    ]


# ============================================================================
# Plotting
# ============================================================================

def generate_plots(rows, output_dir):
    plots_dir = Path(output_dir) / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)

    _, ordered_graphs = build_graph_catalog(rows)
    positions = {g["key"]: i for i, g in enumerate(ordered_graphs)}
    labels = [g["label"] for g in ordered_graphs]
    solvers = sorted({
        solver_label(row) for row in rows if status_ok(row)
    })
    demands = sorted({
        row.get("demand_model")
        for row in rows if row.get("demand_model")
    })

    width = max(14, len(ordered_graphs) * 0.30)
    generated = {
        "congestion": {},
        "runtime": None,
        "failure_worst": None,
        "failure_median": None,
        "failure_critical_25": None,
        "failure_runtime": None,
        "recovery_worst_congestion": None,
        "recovery_degradation": None,
        "recovery_runtime": None,
        "recovery_disconnect_rate": None,
    }

    # Congestion plots
    for demand in demands:
        plt.figure(figsize=(width, 7))
        for solver in solvers:
            points = []
            for row in rows:
                if not status_ok(row):
                    continue
                if solver_label(row) != solver:
                    continue
                if row.get("demand_model") != demand:
                    continue
                value = as_float(row.get("congestion"))
                key = graph_key(row)
                if value is None or value <= 0 or key not in positions:
                    continue
                points.append((positions[key], value))

            points.sort()
            if points:
                plt.plot(
                    [x for x, _ in points],
                    [y for _, y in points],
                    marker="o",
                    markersize=3,
                    linewidth=1.5,
                    label=solver,
                )

        plt.yscale("log")
        plt.xlabel("Graphs ordered by size")
        plt.ylabel("Congestion (log scale)")
        plt.title(f"Congestion by graph size - {demand} demand")
        plt.xticks(range(len(labels)), labels, rotation=90, fontsize=7)
        plt.grid(True, which="both", alpha=0.25)
        plt.legend()
        plt.tight_layout()

        path = plots_dir / f"congestion_{safe_name(demand)}.png"
        plt.savefig(path, dpi=160, bbox_inches="tight")
        plt.close()
        generated["congestion"][demand] = path

    # Solver runtime
    unique_rows = unique_solver_rows_per_graph(rows)
    plt.figure(figsize=(width, 7))
    for solver in solvers:
        points = []
        for row in unique_rows:
            if solver_label(row) != solver:
                continue
            value = as_float(row.get("total_runtime_microseconds"))
            key = graph_key(row)
            if value is None or value <= 0 or key not in positions:
                continue
            points.append((positions[key], value))

        points.sort()
        if points:
            plt.plot(
                [x for x, _ in points],
                [y for _, y in points],
                marker="o",
                markersize=3,
                linewidth=1.5,
                label=solver,
            )

    plt.yscale("log")
    plt.xlabel("Graphs ordered by size")
    plt.ylabel("Runtime [us] (log scale)")
    plt.title("Solver runtime by graph size")
    plt.xticks(range(len(labels)), labels, rotation=90, fontsize=7)
    plt.grid(True, which="both", alpha=0.25)
    plt.legend()
    plt.tight_layout()
    runtime_path = plots_dir / "runtime_by_graph.png"
    plt.savefig(runtime_path, dpi=160, bbox_inches="tight")
    plt.close()
    generated["runtime"] = runtime_path

    if not has_failure_data(rows):
        return generated

    def aggregate_graph_solver(field, reducer):
        grouped = defaultdict(list)
        for row in rows:
            if not status_ok(row):
                continue
            value = as_float(row.get(field))
            if value is None:
                continue
            grouped[(graph_key(row), solver_label(row))].append(value)

        result = {}
        for key, values in grouped.items():
            result[key] = reducer(values)
        return result

    failure_specs = [
        (
            "failure_worst",
            "failure_maximum_lost_traffic_fraction",
            max,
            "Worst single-link traffic exposure by graph",
            "Worst traffic exposure",
            False,
        ),
        (
            "failure_median",
            "failure_median_lost_traffic_fraction",
            statistics.median,
            "Median link exposure by graph",
            "Median traffic exposure",
            False,
        ),
        (
            "failure_critical_25",
            "failure_critical_links_25_percent",
            statistics.median,
            "Critical links (>=25% traffic exposure) by graph",
            "Critical links (median across demands)",
            False,
        ),
        (
            "failure_runtime",
            "failure_analysis_runtime_microseconds",
            statistics.median,
            "Failure-analysis runtime by graph",
            "Runtime [us] (log scale)",
            True,
        ),
    ]

    for output_key, field, reducer, title, ylabel, log_scale in failure_specs:
        values = aggregate_graph_solver(field, reducer)
        plt.figure(figsize=(width, 7))

        for solver in solvers:
            points = []
            for graph in ordered_graphs:
                key = (graph["key"], solver)
                if key not in values:
                    continue
                value = values[key]
                if log_scale and value <= 0:
                    continue
                points.append((positions[graph["key"]], value))

            if points:
                plt.plot(
                    [x for x, _ in points],
                    [y for _, y in points],
                    marker="o",
                    markersize=3,
                    linewidth=1.5,
                    label=solver,
                )

        if log_scale:
            plt.yscale("log")

        plt.xlabel("Graphs ordered by size")
        plt.ylabel(ylabel)
        plt.title(title)
        plt.xticks(range(len(labels)), labels, rotation=90, fontsize=7)
        plt.grid(True, which="both", alpha=0.25)
        plt.legend()
        plt.tight_layout()

        path = plots_dir / f"{output_key}.png"
        plt.savefig(path, dpi=160, bbox_inches="tight")
        plt.close()
        generated[output_key] = path


    if has_recovery_data(rows):
        recovery_specs = [
            (
                "recovery_worst_congestion",
                "recovery_maximum_post_failure_congestion",
                max,
                "Worst survivable post-failure congestion by graph",
                "Post-failure congestion (log scale)",
                True,
            ),
            (
                "recovery_degradation",
                "recovery_maximum_congestion_increase_factor",
                max,
                "Worst survivable congestion degradation by graph",
                "Congestion increase factor",
                False,
            ),
            (
                "recovery_runtime",
                "recovery_average_recomputation_runtime_microseconds",
                statistics.median,
                "N-1 recovery runtime by graph",
                "Average recomputation runtime [us] (log scale)",
                True,
            ),
        ]

        for (
                output_key,
                field,
                reducer,
                title,
                ylabel,
                log_scale,
        ) in recovery_specs:
            values = aggregate_graph_solver(
                field,
                reducer,
            )

            plt.figure(
                figsize=(
                    width,
                    7,
                )
            )

            for solver in solvers:
                points = []

                for graph in ordered_graphs:
                    key = (
                        graph["key"],
                        solver,
                    )

                    if key not in values:
                        continue

                    value = values[key]

                    if (
                            log_scale
                            and value <= 0
                    ):
                        continue

                    points.append(
                        (
                            positions[
                                graph["key"]
                            ],
                            value,
                        )
                    )

                if points:
                    plt.plot(
                        [
                            x
                            for x, _ in points
                        ],
                        [
                            y
                            for _, y in points
                        ],
                        marker="o",
                        markersize=3,
                        linewidth=1.5,
                        label=solver,
                    )

            if log_scale:
                plt.yscale(
                    "log"
                )

            plt.xlabel(
                "Graphs ordered by size"
            )
            plt.ylabel(
                ylabel
            )
            plt.title(
                title
            )
            plt.xticks(
                range(
                    len(labels)
                ),
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
                    /
                    f"{output_key}.png"
            )

            plt.savefig(
                path,
                dpi=160,
                bbox_inches="tight",
            )
            plt.close()

            generated[
                output_key
            ] = path

        disconnect_grouped = defaultdict(list)

        for row in rows:
            if not has_valid_failure_recovery(row):
                continue

            tested = as_float(
                row.get(
                    "recovery_tested_links"
                )
            )
            disconnected = as_float(
                row.get(
                    "recovery_disconnected_failures"
                )
            )

            if (
                    tested is None
                    or tested <= 0
                    or disconnected is None
            ):
                continue

            disconnect_grouped[
                (
                    graph_key(row),
                    solver_label(row),
                )
            ].append(
                disconnected
                /
                tested
            )

        disconnect_values = {
            key:
                statistics.median(values)
            for key, values
            in disconnect_grouped.items()
            if values
        }

        plt.figure(
            figsize=(
                width,
                7,
            )
        )

        for solver in solvers:
            points = []

            for graph in ordered_graphs:
                key = (
                    graph["key"],
                    solver,
                )

                if key not in disconnect_values:
                    continue

                points.append(
                    (
                        positions[
                            graph["key"]
                        ],
                        disconnect_values[
                            key
                        ],
                    )
                )

            if points:
                plt.plot(
                    [
                        x
                        for x, _ in points
                    ],
                    [
                        y
                        for _, y in points
                    ],
                    marker="o",
                    markersize=3,
                    linewidth=1.5,
                    label=solver,
                )

        plt.xlabel(
            "Graphs ordered by size"
        )
        plt.ylabel(
            "Disconnecting link failures / tested links"
        )
        plt.title(
            "Topology disconnection rate under N-1 failures"
        )
        plt.xticks(
            range(
                len(labels)
            ),
            labels,
            rotation=90,
            fontsize=7,
        )
        plt.grid(
            True,
            alpha=0.25,
        )
        plt.legend()
        plt.tight_layout()

        path = (
                plots_dir
                /
                "recovery_disconnect_rate.png"
        )

        plt.savefig(
            path,
            dpi=160,
            bbox_inches="tight",
        )
        plt.close()

        generated[
            "recovery_disconnect_rate"
        ] = path

    return generated


# ============================================================================
# Key findings
# ============================================================================

def generate_key_findings(rows):
    solvers = sorted({
        solver_label(row) for row in rows if status_ok(row)
    })

    congestion = calculate_congestion_summary(rows)
    graph_summary = calculate_graph_summary(rows)

    findings = []

    if solvers:
        congestion_leader = max(
            solvers,
            key=lambda s: congestion[s]["wins"],
        )
        findings.append(
            f"**{congestion_leader}** achieved the lowest congestion most "
            f"frequently, winning or tying **{congestion[congestion_leader]['wins']} "
            f"graph-demand scenarios**."
        )

        runtime_leader = max(
            solvers,
            key=lambda s: graph_summary["runtime_wins"][s],
        )
        median_runtime = median_or_none(
            graph_summary["runtime_values"][runtime_leader]
        )
        findings.append(
            f"**{runtime_leader}** showed the strongest solver-runtime performance, "
            f"finishing fastest on **{graph_summary['runtime_wins'][runtime_leader]} "
            f"graphs** with a median runtime of **{format_runtime_us(median_runtime)}**."
        )

    if has_failure_data(rows):
        failure = calculate_failure_summary(rows)
        failure_solvers = [s for s in solvers if failure[s]["scenarios"]]

        if failure_solvers:
            failure_winner = max(
                failure_solvers,
                key=lambda s: failure[s]["wins"],
            )
            findings.append(
                f"In static single-link exposure analysis, **{failure_winner}** "
                f"most frequently minimized the worst-link traffic loss, winning "
                f"or tying **{failure[failure_winner]['wins']} scenarios**."
            )

            best_median = min(
                failure_solvers,
                key=lambda s: median_or_none(failure[s]["max_loss"]) or float("inf"),
            )
            findings.append(
                f"**{best_median}** had the lowest median worst-link traffic exposure "
                f"at **{format_percentage(median_or_none(failure[best_median]['max_loss']))}**."
            )


    if has_recovery_data(rows):
        recovery = calculate_recovery_summary(
            rows
        )

        recovery_solvers = [
            solver
            for solver in solvers
            if recovery[solver]["scenarios"]
        ]

        if recovery_solvers:
            quality_leader = max(
                recovery_solvers,
                key=lambda solver:
                recovery[
                    solver
                ]["quality_wins"],
            )

            findings.append(
                f"In Layer-2 N-1 recovery analysis, **{quality_leader}** "
                f"most frequently minimized worst survivable post-failure "
                f"congestion, winning or tying **"
                f"{recovery[quality_leader]['quality_wins']} scenarios**."
            )

            speed_leader = max(
                recovery_solvers,
                key=lambda solver:
                recovery[
                    solver
                ]["speed_wins"],
            )

            findings.append(
                f"**{speed_leader}** most frequently delivered the fastest "
                f"post-failure recomputation, winning or tying **"
                f"{recovery[speed_leader]['speed_wins']} scenarios** with "
                f"a median average recovery time of **"
                f"{format_runtime_us(median_or_none(recovery[speed_leader]['avg_runtime']))}**."
            )

    return findings


# ============================================================================
# Report construction
# ============================================================================

def build_report(rows, title, source_path, plots):
    generated = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    _, ordered_graphs = build_graph_catalog(rows)
    solvers = sorted({
        solver_label(row) for row in rows if status_ok(row)
    })

    scenario_groups = {
        (graph_key(row), row.get("demand_model"))
        for row in rows
        if status_ok(row)
    }

    configuration = build_benchmark_configuration(rows)
    congestion_summary = calculate_congestion_summary(rows)
    graph_summary = calculate_graph_summary(rows)

    lines = [
        f"# {title}",
        "",
        f"Generated: `{generated}`",
        "",
        f"Source data: `{source_path}`",
        "",
        "## Benchmark Configuration",
        "",
        markdown_table(
            ["Configuration", "Value"],
            [
                ["Graphs", configuration["graphs"]],
                ["Solvers", configuration["solvers"]],
                ["Demand Models", configuration["demand_models"]],
                ["Random Seed", configuration["seed"]],
                ["Threads", configuration["threads"]],
                ["Result Schema", configuration["schema_version"]],
                ["Layer-2 Failure Recovery", configuration["failure_recovery"]],
            ],
        ),
        "",
        "## Key Findings",
        "",
        (
            f"This benchmark evaluates **{len(solvers)} routing solvers** across "
            f"**{len(ordered_graphs)} unique graph instances** and "
            f"**{len(scenario_groups)} graph-demand scenarios**."
        ),
        "",
    ]

    for finding in generate_key_findings(rows):
        lines.append(f"- {finding}")
    lines.append("")

    # Congestion
    lines += [
        "## Aggregate Congestion Performance",
        "",
        (
            "Congestion is normalized independently for every graph-demand scenario. "
            "A value of `1.0x` means the solver matched the best measured congestion "
            "for that scenario."
        ),
        "",
    ]

    congestion_rows = []
    for solver in solvers:
        data = congestion_summary[solver]
        relative = data["relative"]
        scenarios = data["scenarios"]

        def frac(count):
            return count / scenarios if scenarios else None

        congestion_rows.append([
            solver,
            data["wins"],
            format_number(median_or_none(relative)),
            format_number(geometric_mean(relative)),
            format_percentage(frac(data["within_05"])),
            format_percentage(frac(data["within_10"])),
            format_percentage(frac(data["within_25"])),
        ])

    lines += [
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
            congestion_rows,
        ),
        "",
    ]

    # Solver runtime
    runtime_rows = []
    for solver in solvers:
        runtime_rows.append([
            solver,
            graph_summary["runtime_wins"][solver],
            format_runtime_us(median_or_none(
                graph_summary["runtime_values"][solver]
            )),
            graph_summary["ratio_wins"][solver],
            format_positive_number(
                median_or_none(
                    graph_summary[
                        "ratio_values"
                    ][solver]
                )
            )
        ])

    lines += [
        "## Solver Runtime and Oblivious Ratio",
        "",
        (
            "These metrics are evaluated once per graph because they belong to the "
            "constructed routing scheme and do not depend on the demand model."
        ),
        "",
        markdown_table(
            [
                "Solver",
                "Runtime Wins",
                "Median Solver Runtime",
                "Ratio Wins",
                "Median Oblivious Ratio",
            ],
            runtime_rows,
        ),
        "",
    ]

    # Failure / resilience analysis
    if has_failure_data(rows):
        failure = calculate_failure_summary(rows)

        resilience_rows = []
        for solver in solvers:
            data = failure[solver]
            if not data["scenarios"]:
                continue

            resilience_rows.append([
                solver,
                data["wins"],
                format_percentage(median_or_none(data["max_loss"])),
                format_percentage(median_or_none(data["median_loss"])),
                format_percentage(median_or_none(data["max_affected"])),
                format_number(median_or_none(data["critical_10"])),
                format_number(median_or_none(data["critical_25"])),
                format_number(median_or_none(data["critical_50"])),
                format_runtime_us(median_or_none(data["failure_runtime"])),
                format_runtime_us(median_or_none(data["evaluation_runtime"])),
            ])

        lines += [
            "## Link Failure Exposure Analysis",
            "",
            (
                "**Scope:** this is Layer-1 static single-link exposure analysis. "
                "It measures how much currently routed traffic depends on each physical "
                "link. Actual physical-link removal and fresh post-failure recomputation "
                "are reported separately in the Layer-2 section."
            ),
            "",
            (
                "`Worst-link traffic exposure` is the largest fraction of total routed "
                "traffic that depends on any one tested physical link. Lower is better."
            ),
            "",
            markdown_table(
                [
                    "Solver",
                    "Exposure Wins",
                    "Median Worst-Link Exposure",
                    "Median Link Exposure",
                    "Median Max Affected Demand",
                    "Median Links >=10%",
                    "Median Links >=25%",
                    "Median Links >=50%",
                    "Median Failure Runtime",
                    "Median Demand Eval Runtime",
                ],
                resilience_rows,
            ),
            "",
        ]

        top_rows = []
        for row in top_critical_scenarios(rows, limit=15):
            top_rows.append([
                row.get("graph") or Path(graph_key(row)).stem,
                row.get("demand_model") or "-",
                solver_label(row),
                (
                    f"{row.get('failure_most_critical_source', '?')}"
                    f" -> {row.get('failure_most_critical_target', '?')}"
                ),
                (
                    row.get("failure_most_critical_edge_id")
                    if row.get("failure_most_critical_edge_id") not in (None, "")
                    else "-"
                ),
                format_percentage(
                    as_float(row.get("failure_maximum_lost_traffic_fraction"))
                ),
                format_percentage(
                    as_float(row.get("failure_maximum_affected_demand_fraction"))
                ),
                ])

        lines += [
            "### Highest-Exposure Scenarios",
            "",
            (
                "The following are the highest observed worst-link exposures across "
                "all graph-demand-solver scenarios."
            ),
            "",
            markdown_table(
                [
                    "Graph",
                    "Demand",
                    "Solver",
                    "Critical Link",
                    "Edge ID",
                    "Worst Traffic Exposure",
                    "Demand With Any Exposure",
                ],
                top_rows,
            ),
            "",
        ]

        failure_plot_sections = [
            (
                "### Worst Single-Link Exposure by Graph",
                "failure_worst",
                "For each graph and solver, the plot shows the worst exposure observed across the evaluated demand models.",
            ),
            (
                "### Median Link Exposure by Graph",
                "failure_median",
                "For each graph and solver, the plotted value is the median of the scenario-level median link-exposure values across demand models.",
            ),
            (
                "### Critical Links >=25% Exposure",
                "failure_critical_25",
                "Counts are summarized across demand models for each graph and solver. Lower values indicate fewer highly exposed links.",
            ),
            (
                "### Failure-Analysis Runtime",
                "failure_runtime",
                "This timing is reported separately from solver construction and ordinary demand evaluation.",
            ),
        ]

        for heading, key, description in failure_plot_sections:
            path = plots.get(key)
            if path is None:
                continue
            lines += [
                heading,
                "",
                description,
                "",
                f"![{heading.lstrip('# ').strip()}](plots/{path.name})",
                "",
            ]


    # Layer-2 N-1 recovery analysis
    if has_recovery_data(rows):
        recovery = calculate_recovery_summary(
            rows
        )

        recovery_rows = []

        for solver in solvers:
            data = recovery[solver]

            if not data["scenarios"]:
                continue

            recovery_rows.append([
                solver,
                data["quality_wins"],
                data["speed_wins"],
                format_number(
                    median_or_none(
                        data[
                            "relative_worst_post"
                        ]
                    )
                ),
                format_number(
                    median_or_none(
                        data[
                            "max_post"
                        ]
                    )
                ),
                format_number(
                    median_or_none(
                        data[
                            "avg_post"
                        ]
                    )
                ),
                format_number(
                    median_or_none(
                        data[
                            "max_factor"
                        ]
                    )
                ),
                format_percentage(
                    median_or_none(
                        data[
                            "recoverability_rate"
                        ]
                    )
                ),
                format_percentage(
                    median_or_none(
                        data[
                            "disconnect_rate"
                        ]
                    )
                ),
                format_runtime_us(
                    median_or_none(
                        data[
                            "avg_runtime"
                        ]
                    )
                ),
            ])

        lines += [
            "## Layer-2 N-1 Failure Recovery",
            "",
            (
                "**Scope:** each tested physical link is removed from the topology. "
                "If the remaining topology is connected, the selected solver is "
                "recomputed on the failed topology and the resulting routing is "
                "evaluated on the same demand matrix."
            ),
            "",
            (
                "`Recoverability` is the fraction of survivable N-1 cases for which "
                "solver recomputation succeeded. Disconnecting physical-link failures "
                "are reported separately and are not counted as recomputation failures."
            ),
            "",
            (
                "`Median Worst Recovery vs Best` normalizes the maximum survivable "
                "post-failure congestion independently for every graph-demand scenario. "
                "`1.0x` means the solver matched the best measured recovery quality."
            ),
            "",
            markdown_table(
                [
                    "Solver",
                    "Recovery Quality Wins",
                    "Recovery Speed Wins",
                    "Median Worst Recovery vs Best",
                    "Median Worst Post-Failure Congestion",
                    "Median Avg Post-Failure Congestion",
                    "Median Worst Degradation",
                    "Median Recoverability",
                    "Median Disconnect Rate",
                    "Median Avg Recovery Runtime",
                ],
                recovery_rows,
            ),
            "",
        ]

        top_recovery_rows = []

        for row in top_recovery_scenarios(
                rows,
                limit=15,
        ):
            top_recovery_rows.append([
                row.get("graph")
                or
                Path(
                    graph_key(row)
                ).stem,
                row.get(
                    "demand_model"
                )
                or
                "-",
                solver_label(row),
                format_edge(
                    row.get("recovery_worst_congestion_source"),
                    row.get("recovery_worst_congestion_target"),
                ),
                (
                    row.get(
                        "recovery_worst_congestion_edge_id"
                    )
                    if row.get(
                        "recovery_worst_congestion_edge_id"
                    )
                       not in (
                           None,
                           "",
                       )
                    else "-"
                ),
                format_number(
                    as_float(
                        row.get(
                            "recovery_worst_congestion_post_failure"
                        )
                    )
                ),
                format_number(
                    as_float(
                        row.get(
                            "recovery_worst_congestion_increase_factor"
                        )
                    )
                ),
                format_edge(
                    row.get("recovery_worst_disconnect_source"),
                    row.get("recovery_worst_disconnect_target"),
                ),
                format_percentage(
                    as_float(
                        row.get(
                            "recovery_worst_disconnect_unroutable_demand_fraction"
                        )
                    )
                ),
                format_runtime_us(
                    as_float(
                        row.get(
                            "recovery_slowest_recovery_runtime_microseconds"
                        )
                    )
                ),
                ])

        lines += [
            "### Most Severe Recovery Scenarios",
            "",
            (
                "Scenarios are ordered by their worst survivable congestion "
                "increase factor."
            ),
            "",
            markdown_table(
                [
                    "Graph",
                    "Demand",
                    "Solver",
                    "Worst Survivable Link",
                    "Edge ID",
                    "Post-Failure Congestion",
                    "Degradation",
                    "Worst Disconnect Link",
                    "Unroutable Demand",
                    "Slowest Recovery",
                ],
                top_recovery_rows,
            ),
            "",
        ]

        recovery_plot_sections = [
            (
                "### Worst Survivable Post-Failure Congestion",
                "recovery_worst_congestion",
                (
                    "For each graph and solver, this plot shows the largest "
                    "survivable post-failure congestion observed across demand models."
                ),
            ),
            (
                "### Worst Survivable Congestion Degradation",
                "recovery_degradation",
                (
                    "This plot shows the largest post-failure congestion increase "
                    "factor relative to each solver's own baseline."
                ),
            ),
            (
                "### N-1 Recovery Runtime",
                "recovery_runtime",
                (
                    "Recovery runtime measures fresh post-failure recomputation. "
                    "It is reported separately from baseline solver runtime, demand "
                    "evaluation, and Layer-1 failure-analysis runtime."
                ),
            ),
            (
                "### Topology Disconnection Rate",
                "recovery_disconnect_rate",
                (
                    "This is the fraction of tested physical links whose removal "
                    "disconnects the topology. It is principally a topology property, "
                    "so solvers should agree for the same graph and demand setup."
                ),
            ),
        ]

        for (
                heading,
                key,
                description,
        ) in recovery_plot_sections:
            path = plots.get(
                key
            )

            if path is None:
                continue

            lines += [
                heading,
                "",
                description,
                "",
                (
                    f"![{heading.lstrip('# ').strip()}]"
                    f"(plots/{path.name})"
                ),
                "",
            ]

    # Congestion plots
    lines += [
        "## Congestion by Graph Size",
        "",
        (
            "Graphs are ordered by node count and then edge count. The logarithmic "
            "Y-axis is used because congestion spans several orders of magnitude."
        ),
        "",
    ]

    for demand, path in sorted(plots["congestion"].items()):
        lines += [
            f"### {demand.capitalize()} Demand",
            "",
            f"![Congestion for {demand} demand](plots/{path.name})",
            "",
        ]

    if plots.get("runtime"):
        lines += [
            "## Runtime by Graph Size",
            "",
            (
                "Total solver runtime is shown once per graph and solver. This remains "
                "separate from demand-evaluation runtime and failure-analysis runtime."
            ),
            "",
            f"![Solver runtime by graph size](plots/{plots['runtime'].name})",
            "",
        ]

    # Dataset summary
    dataset = build_dataset_summary(rows)
    lines += [
        "## Benchmark Dataset Summary",
        "",
        markdown_table(
            ["Dataset Property", "Value"],
            [
                ["Graph Instances", dataset["graph_count"]],
                ["Node Range", f"{dataset['min_nodes']} - {dataset['max_nodes']}"],
                ["Edge Range", f"{dataset['min_edges']} - {dataset['max_edges']}"],
            ],
        ),
        "",
        "### Dataset Families",
        "",
        markdown_table(
            ["Dataset Family", "Graphs"],
            [[family, count] for family, count in dataset["families"].items()],
        ),
        "",
        "## Interpretation Notes",
        "",
        "- Lower congestion is better.",
        "- Lower oblivious ratio is better.",
        "- Lower solver runtime is better.",
        "- Lower failure exposure is better.",
        "- Lower post-failure congestion is better.",
        "- Lower recovery runtime is better.",
        (
            "- Solver runtime, demand-evaluation runtime, Layer-1 failure-analysis "
            "runtime, and Layer-2 recomputation runtime are reported separately."
        ),
        (
            "- Layer-1 failure exposure is static dependency analysis. Layer 2 removes "
            "each physical link and performs fresh solver recomputation when the "
            "surviving topology remains connected."
        ),
        (
            "- Disconnecting N-1 failures are reported as topology failures and are "
            "not counted as solver recomputation failures."
        ),
        (
            "- `Demand With Any Exposure` can be close to 100% even when worst traffic "
            "exposure is much lower, because a demand is counted as affected if any "
            "non-zero routed fraction uses the link."
        ),
        (
            "- Full per-instance Layer-1 and Layer-2 recovery measurements are available "
            "in `summary.csv`."
        ),
        "",
    ]

    return "\n".join(lines)


# ============================================================================
# HTML / PDF export
# ============================================================================

def export_report_documents(markdown_path):
    markdown_path = Path(markdown_path)
    html_path = markdown_path.with_suffix(".html")
    pdf_path = markdown_path.with_suffix(".pdf")

    try:
        import markdown
    except Exception as exc:
        print()
        print("[REPORT WARNING] HTML/PDF export skipped because 'markdown' is unavailable.")
        print(f"[REPORT WARNING] {exc}")
        return {"html": None, "pdf": None}

    markdown_text = markdown_path.read_text(encoding="utf-8")
    html_body = markdown.markdown(
        markdown_text,
        extensions=["tables", "fenced_code"],
    )

    css = """
    @page {
        size: A4;
        margin: 16mm;
    }

    body {
        font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Arial, sans-serif;
        font-size: 9.5pt;
        line-height: 1.42;
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
        font-size: 7.5pt;
    }

    th, td {
        border: 1px solid #ccc;
        padding: 4px 6px;
        text-align: left;
        vertical-align: top;
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
        margin: 12px auto 20px auto;
    }

    code {
        font-family: "SFMono-Regular", Consolas, monospace;
        font-size: 0.92em;
        background: #f5f5f5;
        padding: 1px 4px;
    }

    li {
        margin-bottom: 5px;
    }
    """

    complete_html = f"""<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8">
    <title>E-Routing Benchmark Report</title>
    <style>{css}</style>
</head>
<body>
{html_body}
</body>
</html>
"""

    html_path.write_text(complete_html, encoding="utf-8")

    try:
        from weasyprint import HTML
        HTML(
            string=complete_html,
            base_url=str(markdown_path.parent.resolve()),
        ).write_pdf(str(pdf_path))
    except Exception as exc:
        print()
        print("[REPORT WARNING] PDF export skipped.")
        print("[REPORT WARNING] Markdown and HTML reports were still generated.")
        print(f"[REPORT WARNING] {exc}")
        return {"html": html_path, "pdf": None}

    return {"html": html_path, "pdf": pdf_path}


# ============================================================================
# Public API
# ============================================================================

def generate_report_from_summary(
        summary_path,
        report_path=None,
        title="E-Routing Benchmark Report",
):
    summary_path = Path(summary_path)
    if not summary_path.exists():
        raise FileNotFoundError(f"summary.csv does not exist: {summary_path}")

    report_path = (
        Path(report_path)
        if report_path is not None
        else summary_path.parent / "report.md"
    )

    with summary_path.open(newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))

    if not rows:
        raise ValueError("Cannot generate report: summary.csv contains no rows.")

    plots = generate_plots(rows, report_path.parent)
    report = build_report(
        rows=rows,
        title=title,
        source_path=summary_path,
        plots=plots,
    )

    report_path.write_text(report, encoding="utf-8")
    exports = export_report_documents(report_path)

    return {
        "markdown": report_path,
        "html": exports["html"],
        "pdf": exports["pdf"],
    }


# ============================================================================
# CLI
# ============================================================================

def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate E-Routing benchmark, Layer-1 exposure, and Layer-2 N-1 recovery report."
    )
    parser.add_argument("--input", required=True, help="Path to summary.csv.")
    parser.add_argument("--output", default=None, help="Optional report.md output path.")
    parser.add_argument(
        "--title",
        default="E-Routing Benchmark Report",
        help="Report title.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    outputs = generate_report_from_summary(
        summary_path=args.input,
        report_path=args.output,
        title=args.title,
    )

    print()
    print("[REPORT] Generated:")
    print(f"  Markdown: {outputs['markdown']}")
    print(f"  HTML:     {outputs['html'] if outputs['html'] else 'skipped'}")
    print(f"  PDF:      {outputs['pdf'] if outputs['pdf'] else 'skipped'}")


if __name__ == "__main__":
    main()
