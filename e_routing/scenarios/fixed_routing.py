
"""O(m) transformations of validated fixed-routing scenario results."""

from __future__ import annotations

from copy import deepcopy
import math


def transform_fixed_routing(
        baseline: dict,
        *,
        scenario_id: str,
        demand_scale: float,
        capacity_upgrades: tuple[tuple[int, float], ...] = (),
) -> dict:
    if not math.isfinite(demand_scale) or demand_scale <= 0:
        raise ValueError("Invalid demand scale")

    factors = dict(capacity_upgrades)
    if len(factors) != len(capacity_upgrades):
        raise ValueError("Duplicate capacity upgrade")

    if any(
            edge < 0 or not math.isfinite(factor) or factor < 1
            for edge, factor in factors.items()
    ):
        raise ValueError("Invalid capacity upgrade")

    result = deepcopy(baseline)
    links = result["metrics"]["links"]

    if not links:
        raise ValueError("Baseline has no links")

    seen = set()
    maximum = 0.0
    overloaded = []

    for link in links:
        edge_id = link["edge_id"]

        if edge_id in seen:
            raise ValueError("Duplicate baseline edge ID")
        seen.add(edge_id)

        load = link["load"]
        capacity = link["capacity"]

        if (
                not math.isfinite(load)
                or not math.isfinite(capacity)
                or load < 0
                or capacity <= 0
        ):
            raise ValueError("Invalid baseline link metrics")

        factor = factors.get(edge_id, 1.0)

        new_load = load * demand_scale
        new_capacity = capacity * factor
        utilization = new_load / new_capacity

        if (
                not math.isfinite(new_load)
                or not math.isfinite(new_capacity)
                or not math.isfinite(utilization)
                or new_capacity <= 0
        ):
            raise ValueError("Scenario transformation overflow")

        link["load"] = new_load
        link["capacity"] = new_capacity
        link["utilization"] = utilization
        link["overloaded"] = utilization > 1.0

        maximum = max(maximum, utilization)

        if link["overloaded"]:
            overloaded.append(link)

    if not factors.keys() <= seen:
        raise ValueError("Unknown capacity upgrade edge ID")

    result["scenario_id"] = scenario_id
    result["scenario_type"] = "capacity_planning"
    result["execution_mode"] = "fixed_routing_transformation"
    result["demand_scale"] = demand_scale
    result["capacity_upgrades"] = [
        {"edge_id": edge, "capacity_factor": factor}
        for edge, factor in capacity_upgrades
    ]

    result["metrics"]["congestion"] = maximum
    result["metrics"]["overloaded_links"] = overloaded
    result["metrics"]["overloaded_count"] = len(overloaded)
    result["metrics"]["layer1"] = None

    source_evidence = baseline["evidence"]

    result["evidence"] = {
        "graph": source_evidence["graph"],
        "routing_table": source_evidence["routing_table"],
        "derivation": {
            "source_scenario_id": baseline["scenario_id"],
            "source_result_path": source_evidence["result_path"],
            "method": "fixed_routing_load_capacity_transform",
        },
    }

    result["limitations"].append(
        "Derived link metrics use fixed baseline routing and demand; "
        "failure analysis and routing recomputation were not performed."
    )

    return result
