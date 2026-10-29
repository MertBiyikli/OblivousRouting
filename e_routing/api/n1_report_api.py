"""Validate standalone N-1 resilience reports from ResilienceResultJsonWriter.

Keep this independent from the legacy Layer-2 recovery report and the
single-scenario report adapter. Does not execute solvers or access files.
"""

from __future__ import annotations

import math


def _nonnegative(value: object, label: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{label} must be numeric")
    if not math.isfinite(value) or value < 0:
        raise ValueError(f"{label} must be finite and nonnegative")
    return float(value)


def _integer(value: object, label: str) -> int:
    if type(value) is not int or value < 0:
        raise ValueError(f"{label} must be a nonnegative integer")
    return value


def _network(data: object, label: str) -> dict:
    if not isinstance(data, dict):
        raise ValueError(f"{label} must be an object")
    for field in ("max_congestion", "average_congestion", "total_demand_volume"):
        _nonnegative(data.get(field), f"{label}.{field}")
    links = data.get("links")
    if not isinstance(links, list):
        raise ValueError(f"{label}.links must be a list")
    ids = set()
    for link in links:
        if not isinstance(link, dict):
            raise ValueError(f"{label} has invalid link")
        e = _integer(link.get("edge_id"), "edge_id")
        source = _integer(link.get("source"), "source")
        target = _integer(link.get("target"), "target")
        if e in ids or source == target:
            raise ValueError(f"{label} has duplicate edge or self loop")
        ids.add(e)
        capacity = _nonnegative(link.get("capacity"), "capacity")
        load = _nonnegative(link.get("load"), "load")
        utilization = _nonnegative(link.get("utilization"), "utilization")
        if capacity <= 0 or not math.isclose(load / capacity, utilization, rel_tol=1e-7, abs_tol=1e-8):
            raise ValueError(f"{label} inconsistent capacity/load/utilization")
    expected_max = max((x["utilization"] for x in links), default=0.0)
    if not math.isclose(expected_max, data["max_congestion"], rel_tol=1e-7, abs_tol=1e-8):
        raise ValueError(f"{label} max congestion disagrees with links")
    return data


def validate_n1_report(report: object) -> dict:
    if not isinstance(report, dict):
        raise ValueError("N-1 report must be an object")
    if report.get("schema_version") != "1.0" or report.get("analysis_kind") != "n1_resilience":
        raise ValueError("Unsupported N-1 schema")
    fingerprint = report.get("graph_fingerprint")
    if not isinstance(fingerprint, str) or not fingerprint.isdecimal():
        raise ValueError("Invalid graph fingerprint")
    if not isinstance(report.get("solver"), str) or not report["solver"]:
        raise ValueError("Missing solver name")
    baseline = _network(report.get("baseline"), "baseline")
    baseline_edges = {link["edge_id"]: link for link in baseline["links"]}
    outcomes = report.get("outcomes")
    if not isinstance(outcomes, list) or len(outcomes) != len(baseline_edges):
        raise ValueError("N-1 must have one outcome per baseline physical link")
    counts = {"success": 0, "disconnected": 0, "execution_failed": 0}
    seen = set()
    for entry in outcomes:
        if not isinstance(entry, dict):
            raise ValueError("Malformed N-1 outcome")
        eid = _integer(entry.get("edge_id"), "outcome edge_id")
        source = _integer(entry.get("source"), "source")
        target = _integer(entry.get("target"), "target")
        if eid in seen or eid not in baseline_edges:
            raise ValueError("Duplicate or unrecognized failed physical link")
        seen.add(eid)
        base_edge = baseline_edges[eid]
        if (source, target) != (base_edge["source"], base_edge["target"]):
            raise ValueError("Failure endpoints disagree with baseline")
        status = entry.get("status")
        if status not in counts:
            raise ValueError("Invalid N-1 failure status")
        counts[status] += 1
        modified, comparison = entry.get("reoptimized"), entry.get("comparison")
        if status != "success":
            if modified is not None or comparison is not None:
                raise ValueError("Unsuccessful outcome must have null measurements")
            continue
        modified = _network(modified, f"outcome {eid}")
        if {link["edge_id"] for link in modified["links"]} != set(baseline_edges) - {eid}:
            raise ValueError("Modified topology must remove exactly the failed physical link")
        for link in modified["links"]:
            ref = baseline_edges[link["edge_id"]]
            if (link["source"], link["target"]) != (ref["source"], ref["target"]):
                raise ValueError("Surviving physical-link endpoints changed")
        if not math.isclose(modified["total_demand_volume"], baseline["total_demand_volume"], rel_tol=1e-8, abs_tol=1e-8):
            raise ValueError("Demand volume changed across failure")
        if not isinstance(comparison, dict):
            raise ValueError("Missing congestion comparison")
        before = _nonnegative(comparison.get("congestion_before"), "congestion_before")
        after = _nonnegative(comparison.get("congestion_after"), "congestion_after")
        if not math.isclose(before, baseline["max_congestion"], rel_tol=1e-8, abs_tol=1e-8) or not math.isclose(after, modified["max_congestion"], rel_tol=1e-8, abs_tol=1e-8):
            raise ValueError("Comparison and evaluated network disagree")
    summary = report.get("summary")
    if not isinstance(summary, dict) or _integer(summary.get("tested_links"), "tested_links") != len(outcomes):
        raise ValueError("Invalid N-1 summary")
    if any(_integer(summary.get("successful" if key == "success" else key), key) != count for key, count in counts.items()):
        raise ValueError("N-1 summary disagrees with outcomes")
    return report
