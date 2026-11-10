
"""Evidence-driven, single-experiment bottleneck investigation."""

from __future__ import annotations

import json
import math
import re
from pathlib import Path

from e_routing.api.config import load_network_registry
from e_routing.api.scenario_report_api import build_analysis_payload
from e_routing.api.workspace_models import WorkspaceSnapshot
from e_routing.assistant.experiment_executor import (
    run_capacity_experiment,
)
from e_routing.assistant.experiment_models import (
    CapacityExperimentPlan,
)

_ID = re.compile(r"^[A-Za-z0-9_-]{1,64}$")


def select_worst_baseline_link(report: dict) -> dict:
    """Select the highest-utilization physical link.

    The caller must pass a validated detailed-analysis report.
    Ties are resolved using the smallest physical edge ID.
    """
    if report.get("scenario_kind") != "capacity_upgrade":
        raise ValueError("Capacity analysis evidence required")

    solvers = report.get("solver_results")

    if not isinstance(solvers, list) or len(solvers) != 1:
        raise ValueError("Expected one routing-solver result")

    baseline = solvers[0].get("baseline")

    if not isinstance(baseline, dict):
        raise ValueError("Missing baseline network state")

    links = baseline.get("links")

    if not isinstance(links, list) or not links:
        raise ValueError("No baseline links available")

    for link in links:
        if not isinstance(link, dict):
            raise ValueError("Invalid baseline physical link")

        edge_id = link.get("edge_id")
        utilization = link.get("utilization")

        if type(edge_id) is not int or edge_id < 0:
            raise ValueError("Invalid physical edge ID")

        if (
                type(utilization) not in (int, float)
                or not math.isfinite(utilization)
                or utilization < 0
        ):
            raise ValueError("Invalid baseline utilization")

    selected = min(
        links,
        key=lambda link: (
            -link["utilization"],
            link["edge_id"],
        ),
    )

    return {
        "edge_id": selected["edge_id"],
        "source": selected["source"],
        "target": selected["target"],
        "load": selected["load"],
        "capacity": selected["capacity"],
        "utilization": selected["utilization"],
    }


def load_bottleneck_evidence(
        *,
        network_id: str,
        source_analysis_id: str,
        workspace: WorkspaceSnapshot,
        reports_root: Path | str = "results/scenarios",
) -> dict:
    """Load and verify the source analysis before selecting a link."""

    if not _ID.fullmatch(network_id):
        raise ValueError("Invalid network ID")

    if not _ID.fullmatch(source_analysis_id):
        raise ValueError("Invalid source analysis ID")

    registry = load_network_registry()

    network = registry.get(network_id)
    if network is None:
        raise ValueError("Unknown registered network")

    network.verify()

    folder = (
            Path(reports_root).resolve()
            / "detailed"
            / network_id
            / source_analysis_id
    )

    request_path = folder / "request.json"
    report_path = folder / "report.json"

    if not request_path.is_file() or not report_path.is_file():
        raise ValueError("Completed source analysis not found")

    metadata = json.loads(
        request_path.read_text(encoding="utf-8")
    )
    raw_report = json.loads(
        report_path.read_text(encoding="utf-8")
    )

    if not isinstance(metadata, dict):
        raise ValueError("Invalid source request metadata")

    if not isinstance(raw_report, dict):
        raise ValueError("Invalid source analysis report")

    if (
            metadata.get("network_id") != network_id
            or metadata.get("scenario_id") != source_analysis_id
            or metadata.get("scenario_kind") != "capacity_upgrade"
    ):
        raise ValueError("Source analysis identity mismatch")

    if raw_report.get("scenario_kind") != "capacity_upgrade":
        raise ValueError("Source analysis kind mismatch")

    if raw_report.get("status", "completed") != "completed":
        raise ValueError("Source analysis is not completed")

    # Validate all relevant network measurements using the
    # existing trusted detailed-report adapter.
    report = build_analysis_payload(raw_report)

    # The source baseline must correspond to the exact workspace
    # against which the new experiment will be evaluated.
    saved_workspace = WorkspaceSnapshot.model_validate(
        metadata.get("workspace")
    )

    if (
            saved_workspace.model_dump(mode="json")
            != workspace.model_dump(mode="json")
    ):
        raise ValueError(
            "Source analysis workspace differs from the active workspace"
        )

    demand_model = metadata.get("demand_model")

    if demand_model not in (
            "uniform", "gravity", "bimodal", "gaussian"
    ):
        raise ValueError("Invalid source demand model")

    candidate = select_worst_baseline_link(report)

    return {
        "source_analysis_id": source_analysis_id,
        "network_id": network_id,
        "demand_model": demand_model,
        "ranking": "highest_baseline_utilization",
        "candidate": candidate,
        "evaluated_links": len(
            report["solver_results"][0]["baseline"]["links"]
        ),
    }


def run_bottleneck_investigation(
        *,
        network_id: str,
        source_analysis_id: str,
        workspace: WorkspaceSnapshot,
        capacity_factor: float = 1.5,
        reports_root: Path | str = "results/scenarios",
) -> dict:
    """Select one bottleneck and execute one isolated experiment."""

    # Reuse existing executor limitations. Do not silently remove
    # pre-existing upgrades from the comparison workspace.
    if workspace.capacity_upgrades:
        raise ValueError(
            "Bottleneck experiments with existing workspace "
            "capacity upgrades are not yet supported"
        )

    evidence = load_bottleneck_evidence(
        network_id=network_id,
        source_analysis_id=source_analysis_id,
        workspace=workspace,
        reports_root=reports_root,
    )

    candidate = evidence["candidate"]

    plan = CapacityExperimentPlan(
        edge_id=candidate["edge_id"],
        capacity_factor=capacity_factor,
        demand_model=evidence["demand_model"],
    )

    result = run_capacity_experiment(
        network_id=network_id,
        plan=plan,
        workspace=workspace.model_copy(deep=True),
        reports_root=reports_root,
    )

    return {
        "status": "completed",
        "network_id": network_id,
        "source_evidence": evidence,
        "experiment": result.model_dump(mode="json"),
    }
