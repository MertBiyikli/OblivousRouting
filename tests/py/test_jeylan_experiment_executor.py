
"""Deterministic tests for Jeylan capacity experiments."""

from __future__ import annotations

import pytest
from pydantic import ValidationError

import e_routing.assistant.experiment_executor as executor
from e_routing.api.workspace_models import WorkspaceSnapshot
from e_routing.assistant.experiment_models import (
    CapacityExperimentPlan,
)


def test_capacity_experiment_executes_and_preserves_workspace(
        monkeypatch,
        tmp_path,
):
    workspace = WorkspaceSnapshot(demand_scale=1.2)
    original_workspace = workspace.model_dump(mode="json")

    captured = {}

    def fake_execute(*, network_id, request, reports_root):
        captured["network_id"] = network_id
        captured["request"] = request
        captured["reports_root"] = reports_root

        return {
            "network_id": network_id,
            "scenario_id": request.scenario_id,
            "status": "completed",
            "report_available": True,
            "analysis": {
                "scenario_kind": "capacity_upgrade",
                "solver_results": [{
                    "baseline": {"max_congestion": 100.0},
                    "fixed_routing": {"max_congestion": 85.0},
                    "reoptimized": {"max_congestion": 70.0},
                }],
            },
        }

    monkeypatch.setattr(
        executor,
        "execute_detailed_analysis",
        fake_execute,
    )

    result = executor.run_capacity_experiment(
        network_id="test-network",
        plan=CapacityExperimentPlan(
            edge_id=78,
            capacity_factor=1.5,
        ),
        workspace=workspace,
        reports_root=tmp_path,
    )

    request = captured["request"]

    assert captured["network_id"] == "test-network"
    assert captured["reports_root"] == tmp_path

    assert request.scenario_id.startswith("jeylan-")
    assert request.scenario_kind == "capacity_upgrade"
    assert request.link_failures == []

    assert len(request.capacity_upgrades) == 1
    assert request.capacity_upgrades[0].edge_id == 78
    assert request.capacity_upgrades[0].capacity_factor == 1.5

    assert request.workspace.demand_scale == 1.2

    assert result.analysis_id == request.scenario_id
    assert result.evidence_source == "persisted_analysis_report"

    assert result.baseline_max_congestion == 100.0
    assert result.fixed_routing_max_congestion == 85.0
    assert result.reoptimized_max_congestion == 70.0

    assert result.fixed_routing_change == pytest.approx(-0.15)
    assert result.reoptimized_change == pytest.approx(-0.30)

    # No mutation of the caller-owned workspace.
    assert workspace.model_dump(mode="json") == original_workspace
    assert request.workspace is not workspace


@pytest.mark.parametrize(
    "edge_id, factor",
    [
        (-1, 1.5),
        (0, 1.0),
        (78, 101.0),
        (78, float("nan")),
        (78, float("inf")),
    ],
)
def test_invalid_capacity_experiment_is_rejected(
        edge_id,
        factor,
):
    with pytest.raises(ValidationError):
        CapacityExperimentPlan(
            edge_id=edge_id,
            capacity_factor=factor,
        )


def test_existing_capacity_upgrades_are_not_dropped(
        monkeypatch,
):
    def must_not_execute(**kwargs):
        pytest.fail("Execution must not be reached")

    monkeypatch.setattr(
        executor,
        "execute_detailed_analysis",
        must_not_execute,
    )

    workspace = WorkspaceSnapshot(
        capacity_upgrades=[
            {"edge_id": 4, "factor": 1.5},
        ],
    )

    with pytest.raises(
            ValueError,
            match="existing workspace",
    ):
        executor.run_capacity_experiment(
            network_id="test-network",
            plan=CapacityExperimentPlan(
                edge_id=78,
                capacity_factor=1.5,
            ),
            workspace=workspace,
        )


def test_incomplete_experiment_is_rejected(monkeypatch):
    def fake_execute(*, network_id, request, reports_root):
        return {
            "network_id": network_id,
            "scenario_id": request.scenario_id,
            "status": "completed",
            "report_available": False,
        }

    monkeypatch.setattr(
        executor,
        "execute_detailed_analysis",
        fake_execute,
    )

    with pytest.raises(
            RuntimeError,
            match="did not complete",
    ):
        executor.run_capacity_experiment(
            network_id="test-network",
            plan=CapacityExperimentPlan(
                edge_id=78,
                capacity_factor=1.5,
            ),
            workspace=WorkspaceSnapshot(),
        )
