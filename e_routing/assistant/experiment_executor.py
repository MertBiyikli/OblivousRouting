
"""Trusted, simulation-only experiment execution for Jeylan."""

from __future__ import annotations

import math
from pathlib import Path
from uuid import uuid4

from e_routing.api.detailed_analysis_api import (
    DetailedRequest,
    Upgrade,
    execute_detailed_analysis,
)
from e_routing.api.workspace_models import WorkspaceSnapshot
from e_routing.assistant.experiment_models import (
    CapacityExperimentPlan,
    CapacityExperimentResult,
)


def _congestion(value: object, name: str) -> float:
    if type(value) not in (int, float):
        raise ValueError(f"Invalid experiment metric: {name}")

    result = float(value)

    if not math.isfinite(result) or result < 0:
        raise ValueError(f"Invalid experiment metric: {name}")

    return result


def _relative_change(
        baseline: float,
        modified: float,
) -> float | None:
    if baseline <= 0:
        return None

    return (modified - baseline) / baseline


def run_capacity_experiment(
        *,
        network_id: str,
        plan: CapacityExperimentPlan,
        workspace: WorkspaceSnapshot,
        reports_root: Path | str = "results/scenarios",
) -> CapacityExperimentResult:
    """
    Run one isolated capacity-planning experiment.

    Network identity is supplied by the caller's trusted context.
    The active Studio workspace is never mutated.
    """

    # The existing DetailedRequest model does not support
    # a capacity event on top of workspace capacity upgrades.
    # Reject rather than silently discarding those upgrades.
    if workspace.capacity_upgrades:
        raise ValueError(
            "Capacity experiments with existing workspace "
            "capacity upgrades are not yet supported"
        )

    # Server-controlled, collision-resistant experiment ID.
    analysis_id = f"jeylan-{uuid4().hex}"

    request = DetailedRequest(
        scenario_id=analysis_id,
        scenario_kind="capacity_upgrade",
        demand_model=plan.demand_model,
        workspace=workspace.model_copy(deep=True),
        capacity_upgrades=[
            Upgrade(
                edge_id=plan.edge_id,
                capacity_factor=plan.capacity_factor,
            ),
        ],
        link_failures=[],
    )

    # The same service used by Scenario Studio.
    response = execute_detailed_analysis(
        network_id=network_id,
        request=request,
        reports_root=reports_root,
    )

    if not isinstance(response, dict):
        raise RuntimeError("Invalid experiment execution response")

    if (
            response.get("status") != "completed"
            or response.get("scenario_id") != analysis_id
            or response.get("network_id") != network_id
            or response.get("report_available") is not True
    ):
        raise RuntimeError("Experiment did not complete correctly")

    analysis = response.get("analysis")

    if not isinstance(analysis, dict):
        raise RuntimeError("Missing validated experiment report")

    if analysis.get("scenario_kind") != "capacity_upgrade":
        raise RuntimeError("Unexpected experiment analysis type")

    solver_results = analysis.get("solver_results")

    if (
            not isinstance(solver_results, list)
            or len(solver_results) != 1
            or not isinstance(solver_results[0], dict)
    ):
        raise RuntimeError("Invalid solver result structure")

    solver = solver_results[0]

    try:
        baseline = _congestion(
            solver["baseline"]["max_congestion"],
            "baseline.max_congestion",
        )

        fixed = _congestion(
            solver["fixed_routing"]["max_congestion"],
            "fixed_routing.max_congestion",
        )

        reoptimized = _congestion(
            solver["reoptimized"]["max_congestion"],
            "reoptimized.max_congestion",
        )

    except (KeyError, TypeError) as exc:
        raise RuntimeError(
            "Missing experiment congestion measurements"
        ) from exc

    return CapacityExperimentResult(
        network_id=network_id,
        analysis_id=analysis_id,
        edge_id=plan.edge_id,
        capacity_factor=plan.capacity_factor,
        demand_model=plan.demand_model,
        baseline_max_congestion=baseline,
        fixed_routing_max_congestion=fixed,
        reoptimized_max_congestion=reoptimized,
        fixed_routing_change=_relative_change(
            baseline, fixed
        ),
        reoptimized_change=_relative_change(
            baseline, reoptimized
        ),
    )
