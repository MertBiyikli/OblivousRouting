#!/usr/bin/env python3

import argparse
import csv
import math
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description="Compare current routing experiment summary.csv against a baseline."
    )

    parser.add_argument(
        "--baseline",
        required=True,
        help="Path to baseline summary.csv.",
    )

    parser.add_argument(
        "--current",
        required=True,
        help="Path to current summary.csv.",
    )

    parser.add_argument(
        "--congestion-tolerance",
        type=float,
        default=0.05,
        help="Allowed relative congestion increase before regression. Default: 0.05 = 5%%.",
    )

    parser.add_argument(
        "--runtime-tolerance",
        type=float,
        default=0.50,
        help="Allowed relative runtime increase before warning. Default: 0.50 = 50%%.",
    )

    parser.add_argument(
        "--check-runtime",
        action="store_true",
        help="Also compare total_runtime_microseconds and emit warnings.",
    )

    return parser.parse_args()


def normalize(value):
    if value is None:
        return ""

    value = str(value).strip()

    if value.lower() in {"none", "null", "nan"}:
        return ""

    return value


def parse_float(value):
    value = normalize(value)

    if value == "":
        return None

    try:
        x = float(value)
    except ValueError:
        return None

    if not math.isfinite(x):
        return None

    return x


def row_key(row):
    return (
        normalize(row.get("graph")),
        normalize(row.get("solver")),
        normalize(row.get("routing_base")),
        normalize(row.get("demand_model")),
    )


def read_summary(path):
    path = Path(path)

    rows_by_key = {}

    with path.open("r", newline="") as f:
        reader = csv.DictReader(f)

        for row in reader:
            key = row_key(row)

            if key in rows_by_key:
                print(f"[WARN] Duplicate row in {path}: {format_key(key)}")

            rows_by_key[key] = row

    return rows_by_key


def relative_change(baseline, current):
    if baseline == 0:
        if current == 0:
            return 0.0
        return math.inf

    return (current - baseline) / baseline


def format_percent(value):
    if value == math.inf:
        return "+inf"

    return f"{value * 100:+.2f}%"


def format_number(value):
    if value is None:
        return "NA"

    if abs(value) >= 1e6 or (0 < abs(value) < 1e-3):
        return f"{value:.4e}"

    return f"{value:.6g}"


def format_key(key):
    graph, solver, routing_base, demand_model = key

    if routing_base:
        solver_display = f"{solver} [{routing_base}]"
    else:
        solver_display = solver

    return f"{graph} | {solver_display} | {demand_model}"


def compare_metric(
        key,
        baseline_row,
        current_row,
        metric,
        tolerance,
        regression_is_higher=True,
        severity="REGRESSION",
):
    baseline_value = parse_float(baseline_row.get(metric))
    current_value = parse_float(current_row.get(metric))

    if baseline_value is None:
        return {
            "status": "SKIP",
            "severity": "SKIP",
            "key": key,
            "metric": metric,
            "message": f"baseline missing metric {metric}",
        }

    if current_value is None:
        return {
            "status": "MISSING",
            "severity": "REGRESSION",
            "key": key,
            "metric": metric,
            "message": f"current missing metric {metric}",
        }

    change = relative_change(baseline_value, current_value)

    if regression_is_higher:
        is_bad = current_value > baseline_value * (1.0 + tolerance)
    else:
        is_bad = current_value < baseline_value * (1.0 - tolerance)

    if is_bad:
        status = severity
    else:
        status = "PASS"

    return {
        "status": status,
        "severity": severity if is_bad else "PASS",
        "key": key,
        "metric": metric,
        "baseline": baseline_value,
        "current": current_value,
        "change": change,
        "tolerance": tolerance,
    }


def print_result(result):
    status = result["status"]
    key = result["key"]
    metric = result["metric"]

    if status in {"SKIP", "MISSING"}:
        print(f"{status:<10} {format_key(key):<80} {metric:<30} {result['message']}")
        return

    baseline = result["baseline"]
    current = result["current"]
    change = result["change"]

    print(
        f"{status:<10} "
        f"{format_key(key):<80} "
        f"{metric:<30} "
        f"{format_number(baseline):>12} -> {format_number(current):<12} "
        f"{format_percent(change):>10}"
    )


def main():
    args = parse_args()

    baseline = read_summary(args.baseline)
    current = read_summary(args.current)

    results = []

    all_keys = sorted(set(baseline.keys()) | set(current.keys()))

    for key in all_keys:
        baseline_row = baseline.get(key)
        current_row = current.get(key)

        if baseline_row is None:
            results.append({
                "status": "NEW",
                "severity": "PASS",
                "key": key,
                "metric": "-",
                "message": "row exists only in current summary",
            })
            continue

        if current_row is None:
            results.append({
                "status": "MISSING",
                "severity": "REGRESSION",
                "key": key,
                "metric": "-",
                "message": "row missing from current summary",
            })
            continue

        results.append(
            compare_metric(
                key=key,
                baseline_row=baseline_row,
                current_row=current_row,
                metric="congestion",
                tolerance=args.congestion_tolerance,
                regression_is_higher=True,
                severity="REGRESSION",
            )
        )

        if args.check_runtime:
            results.append(
                compare_metric(
                    key=key,
                    baseline_row=baseline_row,
                    current_row=current_row,
                    metric="total_runtime_microseconds",
                    tolerance=args.runtime_tolerance,
                    regression_is_higher=True,
                    severity="WARN",
                )
            )

    print()
    print(f"Baseline: {args.baseline}")
    print(f"Current:  {args.current}")
    print()

    for result in results:
        print_result(result)

    regression_count = sum(
        1 for result in results
        if result.get("severity") == "REGRESSION"
        and result.get("status") in {"REGRESSION", "MISSING"}
    )

    warning_count = sum(
        1 for result in results
        if result.get("severity") == "WARN"
        and result.get("status") == "WARN"
    )

    print()
    print(f"Regressions: {regression_count}")
    print(f"Warnings:    {warning_count}")

    if regression_count > 0:
        print()
        print("Result: REGRESSION")
        raise SystemExit(1)

    if warning_count > 0:
        print()
        print("Result: PASS WITH WARNINGS")
        return

    print()
    print("Result: PASS")


if __name__ == "__main__":
    main()