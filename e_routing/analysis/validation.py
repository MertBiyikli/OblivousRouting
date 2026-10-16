#!/usr/bin/env python3
"""Validate canonical E-Routing JSON against schema + Layer-2 invariants."""
import argparse
import json
import math
from pathlib import Path

from jsonschema import Draft202012Validator, FormatChecker


def near(a, b):
    return math.isclose(a, b, rel_tol=2e-4, abs_tol=2e-4)


def validate_semantics(data):
    errors = []
    n1 = data["n1_recovery"]
    events = n1["events"]
    summary = n1.get("summary")
    tested = n1["tested_link_count"]

    if n1["status"] in ("complete", "completed"):
        if summary is None:
            errors.append("complete N-1 analysis must contain summary")
        if tested is None or tested != len(events):
            errors.append("complete N-1 tested_link_count must equal event count")

    if summary is not None:
        if summary["tested_links"] != tested:
            errors.append("summary.tested_links does not match tested_link_count")
        if len(events) == tested:
            outcomes = {
                "disconnected_failures": sum(e["graph_disconnected"] for e in events),
                "successful_recomputations": sum(e["recomputation_succeeded"] for e in events),
                "failed_recomputations": sum(e["recomputation_attempted"] and not e["recomputation_succeeded"] for e in events),
            }
            for key, expected in outcomes.items():
                if summary[key] != expected:
                    errors.append(f"summary.{key}: expected {expected}, got {summary[key]}")

            measures = (
                ("maximum_post_failure_congestion", "post_failure_congestion"),
                ("maximum_congestion_increase_factor", "congestion_increase_factor"),
                ("maximum_recomputation_runtime_microseconds", "recomputation_runtime_microseconds"),
            )
            for aggregate, event_field in measures:
                available = [e[event_field] for e in events if e[event_field] is not None]
                actual = summary.get(aggregate)
                if available and actual is not None and not near(max(available), actual):
                    errors.append(f"summary.{aggregate} disagrees with event maximum")

    ids = [e["failed_edge_id"] for e in events]
    if len(ids) != len(set(ids)):
        errors.append("duplicate physical-link failure IDs")

    physical = data["topology"]["physical_link_count"]
    if n1["status"] in ("complete", "completed") and physical is not None and tested != physical:
        errors.append("complete N-1 count differs from physical_link_count")

    disconnects = {}
    for e in events:
        edge = e["failed_edge_id"]
        if e["recomputation_succeeded"] and not e["recomputation_attempted"]:
            errors.append(f"edge {edge}: success without attempt")
        if e["graph_disconnected"]:
            disconnects[edge] = e
            if e["recomputation_attempted"]:
                errors.append(f"edge {edge}: recomputation attempted for disconnected graph")
            if e["post_failure_congestion"] is not None:
                errors.append(f"edge {edge}: disconnected failure has post-failure congestion")
        if e["total_demand"] is not None and e["unroutable_demand"] is not None:
            total, lost = e["total_demand"], e["unroutable_demand"]
            if lost > total + 1e-4:
                errors.append(f"edge {edge}: unroutable demand exceeds total")
            if total > 0 and e["unroutable_demand_fraction"] is not None and not near(lost / total, e["unroutable_demand_fraction"]):
                errors.append(f"edge {edge}: unroutable fraction mismatch")
        baseline, post, factor = e["baseline_congestion"], e["post_failure_congestion"], e["congestion_increase_factor"]
        if baseline is not None and baseline > 0 and post is not None and factor is not None and not near(post / baseline, factor):
            errors.append(f"edge {edge}: congestion increase factor mismatch")

    listed = {e["failed_edge_id"]: e for e in data["disconnecting_failures"]}
    if set(listed) != set(disconnects):
        errors.append("disconnecting_failures IDs disagree with N-1 events")
    for edge in set(listed) & set(disconnects):
        a, b = listed[edge], disconnects[edge]
        if (a["source"], a["target"]) != (b["source"], b["target"]):
            errors.append(f"edge {edge}: disconnection endpoints mismatch")
    return errors


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("paths", nargs="+", type=Path, help="JSON files or directories (searched recursively)")
    parser.add_argument("--schema", type=Path, default=Path("schema/analysis-result.schema.json"))
    args = parser.parse_args()
    schema = json.loads(args.schema.read_text())
    Draft202012Validator.check_schema(schema)
    validator = Draft202012Validator(schema, format_checker=FormatChecker())
    files = sorted(set(p for path in args.paths for p in (path.rglob("*.json") if path.is_dir() else [path])))
    if not files:
        parser.error("No JSON files found")
    failures = 0
    for path in files:
        try:
            data = json.loads(path.read_text(), parse_constant=lambda v: (_ for _ in ()).throw(ValueError(f"Invalid JSON constant: {v}")))
            schema_errors = [f"schema at {'/'.join(map(str, err.absolute_path)) or '/'}: {err.message}" for err in validator.iter_errors(data)]
            semantic_errors = validate_semantics(data) if not schema_errors else []
            errors = schema_errors + semantic_errors
        except (ValueError, OSError) as exc:
            errors = [str(exc)]
        if errors:
            failures += 1
            print(f"FAIL {path}")
            for error in errors:
                print(f"  - {error}")
        else:
            print(f"PASS {path}")
    print(f"Validated {len(files)} files: {len(files)-failures} passed, {failures} failed")
    raise SystemExit(1 if failures else 0)


if __name__ == "__main__":
    main()