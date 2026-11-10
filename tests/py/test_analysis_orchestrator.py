
"""Tests for Jeylan's unified network-analysis orchestrator."""

from __future__ import annotations

import asyncio

import pytest
from pydantic import TypeAdapter

import e_routing.assistant.analysis_orchestrator as orchestrator
from e_routing.api.workspace_models import WorkspaceSnapshot
from e_routing.assistant.analysis_plan import AnalysisOperation


ADAPTER = TypeAdapter(AnalysisOperation)


def plan(data):
    return ADAPTER.validate_python(data)


def run(operation, workspace, tmp_path):
    return asyncio.run(
        orchestrator.run_analysis_plan(
            network_id="test-network",
            plan=plan(operation),
            workspace=workspace,
            reports_root=tmp_path,
        )
    )


def test_capacity_dispatch_and_workspace_isolation(
        monkeypatch,
        tmp_path,
):
    calls = []

    def fake_detailed(network_id, request, root):
        calls.append((network_id, request, root))
        return {
            "network_id": network_id,
            "scenario_id": request.scenario_id,
            "status": "completed",
            "analysis": {
                "scenario_kind": "capacity_upgrade",
                "solver_results": [{
                    "baseline": {"max_congestion": 10.0},
                    "fixed_routing": {"max_congestion": 7.0},
                    "reoptimized": {"max_congestion": 8.0},
                }],
            },
        }

    monkeypatch.setattr(
        orchestrator,
        "execute_detailed_analysis",
        fake_detailed,
    )

    workspace = WorkspaceSnapshot(demand_scale=1.2)
    original = workspace.model_dump(mode="json")

    response = run(
        {
            "kind": "capacity_upgrade",
            "upgrades": [
                {"edge_id": 78, "capacity_factor": 1.5}
            ],
        },
        workspace,
        tmp_path,
    )

    assert len(calls) == 1
    request = calls[0][1]

    assert request.scenario_kind == "capacity_upgrade"
    assert request.capacity_upgrades[0].edge_id == 78
    assert request.workspace.demand_scale == 1.2

    assert response["summary"][
               "fixed_routing_max_congestion"
           ] == 7.0

    assert response["analysis_id"].startswith("jeylan-")
    assert response["workspace_modified"] is False
    assert workspace.model_dump(mode="json") == original


def test_link_failure_preserves_existing_upgrades(
        monkeypatch,
        tmp_path,
):
    captured = []

    def fake_detailed(network_id, request, root):
        captured.append(request)
        return {
            "network_id": network_id,
            "scenario_id": request.scenario_id,
            "status": "completed",
            "analysis": {
                "scenario_kind": "link_failure",
                "solver_results": [{
                    "baseline": {"max_congestion": 10.0},
                    "reoptimized": {"max_congestion": 12.0},
                }],
            },
        }

    monkeypatch.setattr(
        orchestrator,
        "execute_detailed_analysis",
        fake_detailed,
    )

    workspace = WorkspaceSnapshot(
        demand_scale=1.2,
        capacity_upgrades=[
            {"edge_id": 4, "factor": 1.5}
        ],
    )

    response = run(
        {
            "kind": "link_failure",
            "failed_edge_id": 14,
        },
        workspace,
        tmp_path,
    )

    request = captured[0]

    assert request.link_failures == [14]
    assert request.workspace.capacity_upgrades[0].edge_id == 4
    assert request.workspace.demand_scale == 1.2
    assert response["summary"][
               "reoptimized_max_congestion"
           ] == 12.0


def test_n1_dispatch_preserves_workspace(
        monkeypatch,
        tmp_path,
):
    calls = []

    async def fake_n1(**kwargs):
        calls.append(kwargs)
        request = kwargs["request"]
        return {
            "network_id": kwargs["network_id"],
            "analysis_id": request.analysis_id,
            "status": "completed",
            "execution_mode": request.execution_mode,
            "summary": {
                "tested_links": 9,
                "successful": 8,
                "disconnected": 1,
                "execution_failed": 0,
            },
        }

    monkeypatch.setattr(
        orchestrator,
        "execute_n1_analysis",
        fake_n1,
    )

    workspace = WorkspaceSnapshot(
        demand_scale=1.3,
        capacity_upgrades=[
            {"edge_id": 4, "factor": 1.5}
        ],
    )

    original = workspace.model_dump(mode="json")

    response = run(
        {
            "kind": "n1_resilience",
            "execution_mode": "fast",
        },
        workspace,
        tmp_path,
    )

    assert len(calls) == 1
    request = calls[0]["request"]

    assert request.workspace.demand_scale == 1.3
    assert len(request.workspace.capacity_upgrades) == 1
    assert calls[0]["reports_root"] == tmp_path / "n1"

    assert response["summary"]["disconnected"] == 1
    assert response["summary"]["successful"] == 8
    assert workspace.model_dump(mode="json") == original


def test_traffic_growth_composes_demand_once(
        monkeypatch,
        tmp_path,
):
    calls = []

    async def fake_traffic(**kwargs):
        calls.append(kwargs)
        request = kwargs["request"]
        return {
            "network_id": kwargs["network_id"],
            "analysis_id": request.analysis_id,
            "status": "completed",
            "execution_mode": "fixed_routing",
            "summary": {
                "baseline_max_congestion": 10.0,
                "modified_max_congestion": 13.0,
            },
        }

    monkeypatch.setattr(
        orchestrator,
        "execute_traffic_analysis",
        fake_traffic,
    )

    workspace = WorkspaceSnapshot(
        demand_scale=1.2,
        capacity_upgrades=[
            {"edge_id": 4, "factor": 1.5}
        ],
        commodity_overrides=[
            {"source": 0, "target": 5, "factor": 2.0}
        ],
    )

    original = workspace.model_dump(mode="json")

    response = run(
        {
            "kind": "traffic_growth",
            "traffic_multiplier": 1.3,
        },
        workspace,
        tmp_path,
    )

    assert len(calls) == 1
    request = calls[0]["request"]

    assert request.global_scale == pytest.approx(1.56)
    assert request.capacity_upgrades[0].edge_id == 4
    assert request.commodity_overrides[0].target == 5
    assert calls[0]["reports_root"] == tmp_path / "traffic"

    assert response["analysis_type"] == "traffic_growth"
    assert workspace.model_dump(mode="json") == original


def test_capacity_rejects_existing_workspace_upgrades(
        monkeypatch,
        tmp_path,
):
    def unexpected_execution(*args, **kwargs):
        pytest.fail("Execution must not be invoked")

    monkeypatch.setattr(
        orchestrator,
        "execute_detailed_analysis",
        unexpected_execution,
    )

    workspace = WorkspaceSnapshot(
        capacity_upgrades=[
            {"edge_id": 4, "factor": 1.5}
        ]
    )

    with pytest.raises(
            ValueError,
            match="existing workspace upgrades",
    ):
        run(
            {
                "kind": "capacity_upgrade",
                "upgrades": [
                    {"edge_id": 78, "capacity_factor": 1.5}
                ],
            },
            workspace,
            tmp_path,
        )


def test_traffic_rejects_excessive_effective_scale(
        monkeypatch,
        tmp_path,
):
    async def unexpected_execution(**kwargs):
        pytest.fail("Traffic service must not execute")

    monkeypatch.setattr(
        orchestrator,
        "execute_traffic_analysis",
        unexpected_execution,
    )

    workspace = WorkspaceSnapshot(demand_scale=100.0)

    with pytest.raises(
            ValueError,
            match="Effective traffic demand multiplier",
    ):
        run(
            {
                "kind": "traffic_growth",
                "traffic_multiplier": 20.0,
            },
            workspace,
            tmp_path,
        )


def test_invalid_network_id_rejected_before_execution(tmp_path):
    with pytest.raises(ValueError, match="network ID"):
        asyncio.run(
            orchestrator.run_analysis_plan(
                network_id="../untrusted",
                plan=plan({
                    "kind": "n1_resilience",
                    "execution_mode": "fast",
                }),
                workspace=WorkspaceSnapshot(),
                reports_root=tmp_path,
            )
        )
