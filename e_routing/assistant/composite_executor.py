
"""Deterministic Traffic Growth + N-1 composite investigation."""

from __future__ import annotations

import json
import logging
import math
import re
from pathlib import Path

from fastapi.concurrency import run_in_threadpool

from e_routing.api.n1_report_api import validate_n1_report
from e_routing.api.workspace_models import WorkspaceSnapshot
from e_routing.assistant.analysis_orchestrator import run_analysis_plan
from e_routing.assistant.composite_plans import (
    TrafficN1InvestigationPlan,
    prepare_traffic_n1_investigation,
)
from e_routing.assistant.composite_evidence import (
    compare_n1_capacity,
)
from e_routing.assistant.composite_plans import (
    CompositeAnalysisPlan,
    CompositeTrafficStep,
    CompositeCapacityStep,
    single_analysis_plan,
)


logger = logging.getLogger(__name__)

_NETWORK_ID = re.compile(r"^[A-Za-z0-9_-]{1,64}$")


def _finite_number(value: object) -> float:
    if type(value) not in (int, float):
        raise ValueError("Invalid numerical report measurement")

    number = float(value)
    if not math.isfinite(number):
        raise ValueError("Non-finite report measurement")

    return number


def _load_verified_n1_report(
        *,
        root: Path,
        network_id: str,
        analysis_id: str,
        effective_workspace: WorkspaceSnapshot,
        demand_model: str,
        execution_mode: str | None = None,
) -> dict:
    """Read a persisted N-1 report and its request provenance."""

    if not _NETWORK_ID.fullmatch(analysis_id):
        raise ValueError("Invalid persisted N-1 analysis ID")

    directory = root / "n1" / network_id / analysis_id

    metadata = json.loads(
        (directory / "request.json").read_text("utf-8")
    )
    raw_report = json.loads(
        (directory / "n1.json").read_text("utf-8")
    )

    if not isinstance(metadata, dict):
        raise ValueError("Invalid persisted N-1 request")

    if (
            metadata.get("network_id") != network_id
            or metadata.get("analysis_id") != analysis_id
    ):
        raise ValueError("N-1 request identity mismatch")

    saved_request = metadata.get("request")
    if not isinstance(saved_request, dict):
        raise ValueError("Missing persisted N-1 configuration")

    if (
            saved_request.get("seed") != 42
            or saved_request.get("demand_model") != demand_model
            or saved_request.get("workspace")
            != effective_workspace.model_dump(mode="json")
    ):
        raise ValueError("N-1 workspace or demand mismatch")

    if (execution_mode is not None and saved_request.get("execution_mode") != execution_mode):
            raise ValueError("N-1 execution mode mismatch")

    report = validate_n1_report(raw_report)

    if report.get("analysis_kind") != "n1_resilience":
        raise ValueError("Unexpected persisted N-1 report type")

    return report


def _rank_failure_outcomes(report: dict) -> dict:
    """Keep disconnections separate from numerical congestion."""

    baseline = report["baseline"]
    before = _finite_number(baseline["max_congestion"])

    disconnected = []
    successful = []
    execution_failed = []
    observed_edges = set()

    for outcome in report["outcomes"]:
        edge_id = outcome["edge_id"]

        if edge_id in observed_edges:
            raise ValueError("Duplicate N-1 physical edge ID")

        observed_edges.add(edge_id)

        entry = {
            "edge_id": edge_id,
            "source": outcome["source"],
            "target": outcome["target"],
            "status": outcome["status"],
        }

        status = outcome["status"]

        if status == "disconnected":
            entry["post_failure_max_congestion"] = None
            entry["relative_congestion_change"] = None
            disconnected.append(entry)

        elif status == "execution_failed":
            entry["post_failure_max_congestion"] = None
            entry["relative_congestion_change"] = None
            execution_failed.append(entry)

        elif status == "success":
            modified = outcome.get("reoptimized")
            if not isinstance(modified, dict):
                raise ValueError(
                    "Successful N-1 failure lacks measurements"
                )

            after = _finite_number(
                modified["max_congestion"]
            )

            entry["post_failure_max_congestion"] = after
            entry["relative_congestion_change"] = (
                (after / before - 1.0)
                if before > 0
                else None
            )
            successful.append(entry)

        else:
            raise ValueError("Unexpected N-1 failure status")

    successful.sort(
        key=lambda entry: (
            -entry["post_failure_max_congestion"],
            entry["edge_id"],
        )
    )

    disconnected.sort(key=lambda item: item["edge_id"])
    execution_failed.sort(key=lambda item: item["edge_id"])

    summary = report["summary"]

    if (
            summary["tested_links"] != len(observed_edges)
            or summary["successful"] != len(successful)
            or summary["disconnected"] != len(disconnected)
            or summary["execution_failed"] != len(execution_failed)
    ):
        raise ValueError("N-1 failure counts do not match report")

    return {
        "baseline_max_congestion": before,
        "tested_links": len(observed_edges),
        "successful": len(successful),
        "disconnected": len(disconnected),
        "execution_failed": len(execution_failed),
        "disconnected_failures": disconnected,
        "execution_failed_failures": execution_failed,
        "worst_successful_failures": successful[:10],
        "ranking_basis": (
            "Post-failure maximum congestion, descending; "
            "physical edge ID breaks ties"
        ),
        "complete_ranking": len(successful) <= 10,
    }


async def execute_traffic_n1_investigation(
        *,
        network_id: str,
        plan: TrafficN1InvestigationPlan,
        workspace: WorkspaceSnapshot,
        reports_root: Path | str = "results/scenarios",
) -> dict:
    """
    Execute Traffic Growth and N-1 under consistent demand.

    Only the two underlying analyses are executed. Their persisted
    reports remain the authoritative measurement sources.
    """

    if (
            not isinstance(network_id, str)
            or not _NETWORK_ID.fullmatch(network_id)
    ):
        raise ValueError("Invalid registered network ID")

    if not isinstance(plan, TrafficN1InvestigationPlan):
        raise ValueError("Invalid composite investigation plan")

    # Current Traffic Growth orchestration uses seed 42.
    if plan.seed != 42:
        raise ValueError(
            "Composite investigations currently require seed 42"
        )

    original = WorkspaceSnapshot.model_validate(
        workspace.model_dump(mode="json")
    ).normalized()

    original_snapshot = original.model_dump(mode="json")

    traffic_plan, n1_plan, effective = (
        prepare_traffic_n1_investigation(
            plan=plan,
            workspace=original,
        )
    )

    root = Path(reports_root).resolve()

    # Traffic receives the original workspace: the orchestrator
    # applies the requested multiplier exactly once.
    traffic = await run_analysis_plan(
        network_id=network_id,
        plan=traffic_plan,
        workspace=original,
        reports_root=root,
    )

    if (
            traffic.get("status") != "completed"
            or traffic.get("analysis_type") != "traffic_growth"
            or traffic.get("network_id") != network_id
    ):
        raise RuntimeError(
            "Traffic Growth execution returned inconsistent evidence"
        )

    traffic_id = traffic["analysis_id"]

    common = {
        "network_id": network_id,
        "kind": "traffic_n1",
        "plan": plan.model_dump(mode="json"),
        "original_workspace": original_snapshot,
        "effective_workspace": effective.model_dump(mode="json"),
        "requested_traffic_multiplier": plan.traffic_multiplier,
        "effective_demand_scale": effective.demand_scale,
        "traffic_analysis_id": traffic_id,
        "traffic_summary": traffic["summary"],
        "workspace_modified": False,
        "evidence_source": "persisted_analysis_reports",
    }

    n1_id = None
    try:
        # N-1 receives the effective workspace directly.
        # Do not apply the traffic multiplier a second time.
        n1 = await run_analysis_plan(
            network_id=network_id,
            plan=n1_plan,
            workspace=effective,
            reports_root=root,
        )

        if isinstance(n1, dict):
            n1_id = n1.get("analysis_id")

        if (
                n1.get("status") != "completed"
                or n1.get("analysis_type") != "n1_resilience"
                or n1.get("network_id") != network_id
        ):
            raise ValueError("N-1 execution identity mismatch")

        report = await run_in_threadpool(
            _load_verified_n1_report,
            root=root,
            network_id=network_id,
            analysis_id=n1["analysis_id"],
            effective_workspace=effective,
            demand_model=plan.demand_model,
        )

        ranking = _rank_failure_outcomes(report)

        # Compare actual generated demand volumes, not the maximum
        # congestion: Traffic uses fixed routing, N-1 recomputes it.
        traffic_result = traffic.get("service_result", {})
        traffic_summary = traffic_result.get("summary", {})

        traffic_volume = traffic_summary.get(
            "modified_demand_volume"
        )

        if traffic_volume is None:
            raise ValueError(
                "Traffic Growth demand-volume evidence is missing"
            )

        observed_volume = _finite_number(
            report["baseline"]["total_demand_volume"]
        )
        expected_volume = _finite_number(traffic_volume)

        if not math.isclose(
                observed_volume,
                expected_volume,
                rel_tol=1e-6,
                abs_tol=1e-6,
        ):
            raise ValueError(
                "Traffic and N-1 demand volumes differ"
            )

        return {
            **common,
            "status": "completed",
            "n1_analysis_id": n1_id,
            "n1_summary": n1["summary"],
            "resilience": ranking,
            "n1_resilient": (
                    ranking["disconnected"] == 0
                    and ranking["execution_failed"] == 0
            ),
            "congestion_limit_checked": False,
            "capacity_candidates_tested": 0,
        }

    except Exception:
        # The first experiment is already persisted.
        # Preserve its ID rather than losing successful work.
        logger.exception(
            "Composite investigation incomplete after "
            "Traffic Growth analysis %s",
            traffic_id,
        )

        return {
            **common,
            "status": "partial",
            "n1_analysis_id": (
                n1.get("analysis_id")
                if isinstance(locals().get("n1"), dict)
                else None
            ),
            "n1_summary": None,
            "resilience": None,
            "n1_resilient": None,
            "congestion_limit_checked": False,
            "capacity_candidates_tested": 0,
            "error": (
                "N-1 execution or persisted evidence validation "
                "did not complete successfully"
            ),
        }



# ------------------------------------------------------------
# Generalized sequential composite executor — F3-A
# ------------------------------------------------------------

def _transition_workspace(
        workspace: WorkspaceSnapshot,
        step,
) -> WorkspaceSnapshot:
    """Calculate the next simulated workspace."""

    data = workspace.model_dump(mode="json")

    if isinstance(step, CompositeTrafficStep):
        scale = (
                workspace.demand_scale
                * step.traffic_multiplier
        )

        if not math.isfinite(scale) or not 0 < scale <= 1000:
            raise ValueError(
                "Effective demand scale is out of range"
            )

        data["demand_scale"] = scale

    elif isinstance(step, CompositeCapacityStep):
        if data["capacity_upgrades"]:
            raise ValueError(
                "Capacity comparison requires a workspace "
                "without existing capacity upgrades"
            )

        for upgrade in step.upgrades:
            data["capacity_upgrades"].append({
                "edge_id": upgrade.edge_id,
                "factor": upgrade.capacity_factor,
            })

    return WorkspaceSnapshot.model_validate(
        data
    ).normalized()


def _preflight_composite(
        plan: CompositeAnalysisPlan,
        workspace: WorkspaceSnapshot,
) -> None:
    """Validate the whole workspace transition sequence."""

    simulated = workspace

    for step in plan.steps:
        # Reject unsupported transitions before any execution.
        simulated = _transition_workspace(
            simulated,
            step,
        )



async def _build_composite_comparisons(
        *,
        network_id: str,
        completed_steps: list[dict],
        reports_root: Path,
) -> list[dict]:
    """Compare consecutive N-1 observations with capacity changes.

    Other evaluation steps may occur between N-1 observations,
    but demand settings must remain identical.
    """

    n1_steps = [
        step
        for step in completed_steps
        if step["operation"] == "n1_resilience"
    ]

    comparisons = []

    for before, after in zip(n1_steps, n1_steps[1:]):
        before_ws = WorkspaceSnapshot.model_validate(
            before["input_workspace"]
        ).normalized()

        after_ws = WorkspaceSnapshot.model_validate(
            after["input_workspace"]
        ).normalized()

        # Comparison only applies to a capacity configuration
        # change with otherwise identical demand.
        if (
                before_ws.demand_scale != after_ws.demand_scale
                or before_ws.commodity_overrides
                != after_ws.commodity_overrides
        ):
            continue

        if (
                before_ws.capacity_upgrades
                == after_ws.capacity_upgrades
        ):
            continue

        before_plan = before["plan"]
        after_plan = after["plan"]

        if (
                before_plan["demand_model"]
                != after_plan["demand_model"]
                or before_plan["seed"] != after_plan["seed"]
                or before_plan["execution_mode"]
                != after_plan["execution_mode"]
        ):
            continue

        old_report = await run_in_threadpool(
            _load_verified_n1_report,
            root=reports_root,
            network_id=network_id,
            analysis_id=before["analysis_id"],
            effective_workspace=before_ws,
            demand_model=before_plan["demand_model"],
            execution_mode=before_plan["execution_mode"],
        )

        new_report = await run_in_threadpool(
            _load_verified_n1_report,
            root=reports_root,
            network_id=network_id,
            analysis_id=after["analysis_id"],
            effective_workspace=after_ws,
            demand_model=after_plan["demand_model"],
            execution_mode=after_plan["execution_mode"],
        )

        comparison = compare_n1_capacity(
            before_report=old_report,
            after_report=new_report,
            before_workspace=before_ws,
            after_workspace=after_ws,
            before_analysis_id=before["analysis_id"],
            after_analysis_id=after["analysis_id"],
        )

        comparisons.append(comparison)

    return comparisons


async def execute_composite_analysis(
        *,
        network_id: str,
        plan: CompositeAnalysisPlan,
        workspace: WorkspaceSnapshot,
        reports_root: Path | str = "results/scenarios",
) -> dict:
    """Execute validated operations sequentially.

    Individual analyses are delegated to the trusted
    orchestrator. Only completed steps advance the
    simulated workspace.
    """

    if (
            not isinstance(network_id, str)
            or not _NETWORK_ID.fullmatch(network_id)
    ):
        raise ValueError("Invalid network ID")

    if not isinstance(plan, CompositeAnalysisPlan):
        raise ValueError("Invalid composite plan")

    if not isinstance(workspace, WorkspaceSnapshot):
        raise ValueError("Invalid workspace")

    original = WorkspaceSnapshot.model_validate(
        workspace.model_dump(mode="json")
    ).normalized()

    _preflight_composite(plan, original)

    current = original
    completed = []
    root = Path(reports_root).resolve()

    for index, step in enumerate(plan.steps):
        before = current.model_dump(mode="json")

        try:
            operation_plan = single_analysis_plan(
                step,
                plan,
            )

            result = await run_analysis_plan(
                network_id=network_id,
                plan=operation_plan,
                workspace=current,
                reports_root=root,
            )

            if (
                    result.get("status") != "completed"
                    or result.get("network_id") != network_id
                    or result.get("analysis_type")
                    != step.operation
            ):
                raise RuntimeError(
                    "Inconsistent analysis execution result"
                )

            next_workspace = _transition_workspace(
                current,
                step,
            )

        except Exception:
            logger.exception(
                "Composite step %s failed", index
            )

            if not completed:
                raise

            return {
                "network_id": network_id,
                "kind": "composite",
                "status": "partial",
                "plan": plan.model_dump(mode="json"),
                "original_workspace":
                    original.model_dump(mode="json"),
                "effective_workspace":
                    current.model_dump(mode="json"),
                "completed_steps": completed,
                "failed_step_index": index,
                "failed_operation": step.operation,
                "workspace_modified": False,
                "error": (
                    "Composite execution stopped. "
                    "Completed reports remain persisted; "
                    "inspect reports before retrying."
                ),
            }

        completed.append({
            "step_index": index,
            "operation": step.operation,
            "analysis_id": result["analysis_id"],
            "status": "completed",
            "plan": result["plan"],
            "summary": result["summary"],
            "input_workspace": before,
            "output_workspace":
                next_workspace.model_dump(mode="json"),
            "evidence_source":
                "persisted_analysis_report",
        })

        current = next_workspace


    try:
        comparisons = await _build_composite_comparisons(
            network_id=network_id,
            completed_steps=completed,
            reports_root=root,
        )

    except Exception:
        logger.exception(
            "Composite comparison verification failed"
        )

        return {
            "network_id": network_id,
            "kind": "composite",
            "status": "partial",
            "plan": plan.model_dump(mode="json"),
            "original_workspace":
                original.model_dump(mode="json"),
            "effective_workspace":
                current.model_dump(mode="json"),
            "completed_steps": completed,
            "failed_step_index": None,
            "failed_operation": None,
            "comparisons": [],
            "comparison_status": "failed",
            "workspace_modified": False,
            "error": (
                "All requested analyses completed, but "
                "cross-analysis comparison validation failed. "
                "Persisted analysis reports remain available."
            ),
        }

    return {
        "network_id": network_id,
        "kind": "composite",
        "status": "completed",
        "plan": plan.model_dump(mode="json"),
        "original_workspace":
            original.model_dump(mode="json"),
        "effective_workspace":
            current.model_dump(mode="json"),
        "completed_steps": completed,
        "failed_step_index": None,
        "comparisons": comparisons,
        "comparison_status": "completed",
        "workspace_modified": False,
        "evidence_source":
            "persisted_analysis_reports",
    }

