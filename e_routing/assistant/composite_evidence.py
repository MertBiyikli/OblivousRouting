
"""Deterministic comparisons of verified composite N-1 evidence."""

from __future__ import annotations

import math

from e_routing.api.workspace_models import WorkspaceSnapshot

REL_TOL = 1e-9
ABS_TOL = 1e-9

STATUSES = {"success", "disconnected", "execution_failed"}


def _finite(value: object) -> float:
    if type(value) not in (int, float):
        raise ValueError("Invalid numeric measurement")

    number = float(value)

    if not math.isfinite(number):
        raise ValueError("Non-finite numeric measurement")

    return number


def _relative(before: float, after: float) -> float | None:
    return (after / before - 1.0) if before > 0 else None


def _outcomes(report: dict) -> dict:
    mapped = {}

    for outcome in report["outcomes"]:
        edge_id = outcome["edge_id"]

        if type(edge_id) is not int or edge_id < 0:
            raise ValueError("Invalid physical edge ID")

        if edge_id in mapped:
            raise ValueError("Duplicate physical edge ID")

        if outcome["status"] not in STATUSES:
            raise ValueError("Unknown N-1 status")

        mapped[edge_id] = outcome

    if report["summary"]["tested_links"] != len(mapped):
        raise ValueError("N-1 outcome count mismatch")

    return mapped


def _congestion(outcome: dict) -> float | None:
    if outcome["status"] != "success":
        return None

    measurements = outcome.get("reoptimized")

    if not isinstance(measurements, dict):
        raise ValueError("Missing successful failure measurements")

    return _finite(measurements["max_congestion"])


def compare_n1_capacity(
        *,
        before_report: dict,
        after_report: dict,
        before_workspace: WorkspaceSnapshot,
        after_workspace: WorkspaceSnapshot,
        before_analysis_id: str,
        after_analysis_id: str,
) -> dict:
    """Compare verified N-1 reports with identical demand settings.

    The caller is responsible for validating persisted report
    identity, demand model, seed and execution mode.
    """

    before_ws = before_workspace.normalized()
    after_ws = after_workspace.normalized()

    if (
            before_ws.demand_scale != after_ws.demand_scale
            or before_ws.commodity_overrides
            != after_ws.commodity_overrides
    ):
        raise ValueError("Demand configurations differ")

    if before_ws.capacity_upgrades == after_ws.capacity_upgrades:
        raise ValueError("Capacity configurations are identical")

    old_outcomes = _outcomes(before_report)
    new_outcomes = _outcomes(after_report)

    if old_outcomes.keys() != new_outcomes.keys():
        raise ValueError("Physical edge sets differ")

    old_volume = _finite(
        before_report["baseline"]["total_demand_volume"]
    )
    new_volume = _finite(
        after_report["baseline"]["total_demand_volume"]
    )

    if not math.isclose(
            old_volume, new_volume, rel_tol=1e-6, abs_tol=1e-6
    ):
        raise ValueError("Demand volumes differ")

    old_baseline = _finite(
        before_report["baseline"]["max_congestion"]
    )
    new_baseline = _finite(
        after_report["baseline"]["max_congestion"]
    )

    counts = {
        "improved": 0,
        "regressed": 0,
        "unchanged": 0,
        "became_connected": 0,
        "became_disconnected": 0,
        "execution_status_changed": 0,
        "execution_failed_both": 0,
        "disconnected_both": 0,
    }

    comparisons = []

    for edge_id in sorted(old_outcomes):
        old = old_outcomes[edge_id]
        new = new_outcomes[edge_id]

        if (
                old["source"] != new["source"]
                or old["target"] != new["target"]
        ):
            raise ValueError("Physical edge endpoints differ")

        old_status = old["status"]
        new_status = new["status"]

        old_value = _congestion(old)
        new_value = _congestion(new)

        change = None

        if "execution_failed" in (old_status, new_status):
            classification = (
                "execution_failed_both"
                if old_status == new_status
                else "execution_status_changed"
            )

        elif old_status == new_status == "disconnected":
            classification = "disconnected_both"

        elif old_status == "disconnected":
            classification = "became_connected"

        elif new_status == "disconnected":
            classification = "became_disconnected"

        else:
            if old_value is None or new_value is None:
                raise ValueError("Missing congestion values")

            change = _relative(old_value, new_value)

            if math.isclose(
                    old_value,
                    new_value,
                    rel_tol=REL_TOL,
                    abs_tol=ABS_TOL,
            ):
                classification = "unchanged"
            elif new_value < old_value:
                classification = "improved"
            else:
                classification = "regressed"

        counts[classification] += 1

        comparisons.append({
            "edge_id": edge_id,
            "source": old["source"],
            "target": old["target"],
            "before_status": old_status,
            "after_status": new_status,
            "before_max_congestion": old_value,
            "after_max_congestion": new_value,
            "relative_congestion_change": change,
            "classification": classification,
        })

    assert sum(counts.values()) == len(comparisons)

    return {
        "comparison_type": "capacity_n1",
        "before_analysis_id": before_analysis_id,
        "after_analysis_id": after_analysis_id,
        "before_capacity_upgrades": [
            item.model_dump(mode="json")
            for item in before_ws.capacity_upgrades
        ],
        "after_capacity_upgrades": [
            item.model_dump(mode="json")
            for item in after_ws.capacity_upgrades
        ],
        "demand_scale": before_ws.demand_scale,
        "total_demand_volume": old_volume,
        "before_baseline_max_congestion": old_baseline,
        "after_baseline_max_congestion": new_baseline,
        "relative_baseline_congestion_change": _relative(
            old_baseline, new_baseline
        ),
        "counts": counts,
        "edge_comparisons": comparisons,
        "evidence_source": "persisted_n1_reports",
    }
