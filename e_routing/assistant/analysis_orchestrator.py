
"""Trusted dispatcher for Jeylan's four network-analysis operations."""

from __future__ import annotations

import math
import re
from pathlib import Path
from uuid import uuid4

from fastapi.concurrency import run_in_threadpool

from e_routing.api.detailed_analysis_api import (
    DetailedRequest,
    Upgrade,
    execute_detailed_analysis,
)
from e_routing.api.n1_service import (
    N1Request,
    execute_n1_analysis,
)
from e_routing.api.traffic_models import (
    TrafficRequest,
)
from e_routing.api.traffic_service import (
    execute_traffic_analysis,
)
from e_routing.api.workspace_models import WorkspaceSnapshot
from e_routing.assistant.analysis_plan import (
    AnalysisOperation,
    CapacityAnalysisPlan,
    LinkFailureAnalysisPlan,
    N1AnalysisPlan,
    TrafficGrowthAnalysisPlan,
)


_NETWORK_ID = re.compile(r"^[A-Za-z0-9_-]{1,64}$")


def _new_analysis_id() -> str:
    """Generate an analysis ID independent of model output."""
    return f"jeylan-{uuid4().hex}"


def _validated_workspace(
        workspace: WorkspaceSnapshot,
) -> WorkspaceSnapshot:
    """Validate and isolate the caller's workspace."""
    if not isinstance(workspace, WorkspaceSnapshot):
        raise ValueError("Invalid network workspace")

    return WorkspaceSnapshot.model_validate(
        workspace.model_dump(mode="json")
    ).normalized()


def _effective_traffic_workspace(
        workspace: WorkspaceSnapshot,
        multiplier: float,
) -> WorkspaceSnapshot:
    """Compose relative traffic growth with existing demand."""
    effective_scale = workspace.demand_scale * multiplier

    if (
            not math.isfinite(effective_scale)
            or effective_scale <= 0
            or effective_scale > 1000
    ):
        raise ValueError(
            "Effective traffic demand multiplier is out of range"
        )

    data = workspace.model_dump(mode="json")
    data["demand_scale"] = effective_scale

    return WorkspaceSnapshot.model_validate(data).normalized()


def _detailed_summary(result: dict) -> dict:
    """Extract measured congestion if a detailed report is present."""
    report = result.get("analysis")

    if not isinstance(report, dict):
        return {}

    solvers = report.get("solver_results")

    if not isinstance(solvers, list) or len(solvers) != 1:
        return {}

    solver = solvers[0]
    baseline = solver.get("baseline")
    fixed = solver.get("fixed_routing")
    reoptimized = solver.get("reoptimized")

    if not isinstance(baseline, dict):
        return {}

    return {
        "baseline_max_congestion": baseline.get(
            "max_congestion"
        ),
        "fixed_routing_max_congestion": (
            fixed.get("max_congestion")
            if isinstance(fixed, dict)
            else None
        ),
        "reoptimized_max_congestion": (
            reoptimized.get("max_congestion")
            if isinstance(reoptimized, dict)
            else None
        ),
        "disconnected": (
                result.get("analysis", {}).get("scenario_kind")
                == "link_failure"
                and reoptimized is None
        ),
    }


async def run_analysis_plan(
        *,
        network_id: str,
        plan: AnalysisOperation,
        workspace: WorkspaceSnapshot,
        reports_root: Path | str,
) -> dict:
    """Execute one validated network analysis.

    The caller is responsible for authorization, concurrency
    limits and the execution budget.

    This function uses only registered networks and
    server-controlled report destinations.
    """

    if not isinstance(network_id, str) or not _NETWORK_ID.fullmatch(
            network_id
    ):
        raise ValueError("Invalid registered network ID")

    if not isinstance(
            plan,
            (
                    CapacityAnalysisPlan,
                    LinkFailureAnalysisPlan,
                    N1AnalysisPlan,
                    TrafficGrowthAnalysisPlan,
            ),
    ):
        raise ValueError("Unsupported analysis plan")

    # Take a validated, independent snapshot before execution.
    original_workspace = _validated_workspace(workspace)

    root = Path(reports_root).resolve()
    analysis_id = _new_analysis_id()

    # ---------------------------------------------------------
    # Capacity Upgrade
    # ---------------------------------------------------------

    if isinstance(plan, CapacityAnalysisPlan):
        if original_workspace.capacity_upgrades:
            raise ValueError(
                "Three-way capacity comparison with existing "
                "workspace upgrades is not supported"
            )

        request = DetailedRequest(
            scenario_id=analysis_id,
            scenario_kind="capacity_upgrade",
            demand_model=plan.demand_model,
            workspace=original_workspace,
            capacity_upgrades=[
                Upgrade(
                    edge_id=item.edge_id,
                    capacity_factor=item.capacity_factor,
                )
                for item in plan.upgrades
            ],
        )

        result = await run_in_threadpool(
            execute_detailed_analysis,
            network_id,
            request,
            root,
        )

        summary = _detailed_summary(result)

    # ---------------------------------------------------------
    # Link Failure
    # ---------------------------------------------------------

    elif isinstance(plan, LinkFailureAnalysisPlan):
        request = DetailedRequest(
            scenario_id=analysis_id,
            scenario_kind="link_failure",
            demand_model=plan.demand_model,
            workspace=original_workspace,
            link_failures=[plan.failed_edge_id],
        )

        result = await run_in_threadpool(
            execute_detailed_analysis,
            network_id,
            request,
            root,
        )

        summary = _detailed_summary(result)

    # ---------------------------------------------------------
    # N-1 Resilience
    # ---------------------------------------------------------

    elif isinstance(plan, N1AnalysisPlan):
        request = N1Request(
            analysis_id=analysis_id,
            demand_model=plan.demand_model,
            seed=plan.seed,
            execution_mode=plan.execution_mode,
            workspace=original_workspace,
        )

        result = await execute_n1_analysis(
            network_id=network_id,
            request=request,
            reports_root=root / "n1",
        )

        summary = result["summary"]

    # ---------------------------------------------------------
    # Traffic Growth
    # ---------------------------------------------------------

    elif isinstance(plan, TrafficGrowthAnalysisPlan):
        effective_workspace = _effective_traffic_workspace(
            original_workspace,
            plan.traffic_multiplier,
        )

        request = TrafficRequest(
            analysis_id=analysis_id,
            demand_model=plan.demand_model,
            global_scale=effective_workspace.demand_scale,
            capacity_upgrades=[
                {
                    "edge_id": item.edge_id,
                    "capacity_factor": item.factor,
                }
                for item in effective_workspace.capacity_upgrades
            ],
            commodity_overrides=[
                {
                    "source": item.source,
                    "target": item.target,
                    "factor": item.factor,
                }
                for item in effective_workspace.commodity_overrides
            ],
        )

        result = await execute_traffic_analysis(
            network_id=network_id,
            request=request,
            reports_root=root / "traffic",
        )

        summary = result["summary"]

    else:
        raise ValueError("Unsupported analysis plan")

    # The backend owns the returned identity. Never take it
    # from the language model or an unverified report field.
    if (
            result.get("network_id") != network_id
            or result.get("status") != "completed"
    ):
        raise RuntimeError(
            "Analysis service returned an inconsistent result"
        )

    returned_id = result.get(
        "scenario_id",
        result.get("analysis_id"),
    )

    if returned_id != analysis_id:
        raise RuntimeError(
            "Analysis service returned an unexpected analysis ID"
        )

    return {
        "network_id": network_id,
        "analysis_id": analysis_id,
        "analysis_type": plan.kind,
        "status": "completed",
        "summary": summary,
        "plan": plan.model_dump(mode="json"),
        "evidence_source": "persisted_analysis_report",
        "workspace_modified": False,
        "service_result": result,
    }
