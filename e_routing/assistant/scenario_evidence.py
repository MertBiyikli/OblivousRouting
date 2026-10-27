
"""Read-only evidence extraction for persisted-routing scenarios."""

from __future__ import annotations
import math
from typing import Any


def summarize_scenario(result: dict[str, Any]) -> dict[str, Any]:
    if result.get("status") != "completed":
        raise ValueError("Scenario must be completed")

    scenario_id = result.get("scenario_id")
    metrics = result.get("metrics")

    if not isinstance(scenario_id, str) or not scenario_id:
        raise ValueError("Scenario ID is missing")
    if not isinstance(metrics, dict):
        raise ValueError("Scenario metrics are missing")

    links = metrics.get("links")
    if not isinstance(links, list):
        raise ValueError("Scenario link evidence is missing")

    def valid_ratio(value: Any) -> bool:
            return type(value) in (int, float) and math.isfinite(value) and value >= 0

    congestion = metrics.get("congestion")
    if not valid_ratio(congestion):
        raise ValueError("Invalid scenario congestion")

    overloaded_count = metrics.get("overloaded_count")
    if type(overloaded_count) is not int or overloaded_count < 0:
        raise ValueError("Invalid overloaded-link count")

    seen_edges = set()
    actual_overloaded = []

    for link in links:
        if not isinstance(link, dict):
            raise ValueError("Invalid link evidence")

        edge_id = link.get("edge_id")
        utilization = link.get("utilization")
        overloaded = link.get("overloaded")

        if type(edge_id) is not int or edge_id < 0 or edge_id in seen_edges:
            raise ValueError("Invalid or duplicate edge ID")

        if not valid_ratio(utilization) or type(overloaded) is not bool:
            raise ValueError("Invalid link utilization evidence")

        seen_edges.add(edge_id)

        if overloaded != (utilization > 1.0):
            raise ValueError("Inconsistent overloaded-link classification")

        if overloaded:
            actual_overloaded.append(link)

    if overloaded_count != len(actual_overloaded):
        raise ValueError("Overloaded-link count mismatch")

    if metrics.get("overloaded_links") != actual_overloaded:
        raise ValueError("Overloaded-link evidence mismatch")

    return {
        "evidence_source": "persisted_routing_scenario",
        "scenario_id": scenario_id,
        "execution_mode": result.get("execution_mode"),
        "scenario_type": result.get("scenario_type"),
        "demand_model": result.get("demand_model"),
        "demand_scale": result.get("demand_scale"),
        "capacity_upgrades": result.get("capacity_upgrades"),
        "congestion_ratio": metrics.get("congestion"),
        "overloaded_count": metrics.get("overloaded_count"),
        "overloaded_links": metrics.get("overloaded_links"),
        "links": links,
        "layer1": metrics.get("layer1"),
        "limitations": result.get("limitations", []),
    }

