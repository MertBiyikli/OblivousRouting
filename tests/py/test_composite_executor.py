
"""Regression tests for composite Traffic Growth + N-1."""

import asyncio
from copy import deepcopy

import pytest

import e_routing.assistant.composite_executor as composite

from e_routing.api.workspace_models import WorkspaceSnapshot
from e_routing.assistant.composite_plans import (
    TrafficN1InvestigationPlan,
)


def sample_n1_report():
    return {
        "baseline": {
            "max_congestion": 200.0,
            "total_demand_volume": 1560.0,
        },
        "summary": {
            "tested_links": 3,
            "successful": 2,
            "disconnected": 1,
            "execution_failed": 0,
        },
        "outcomes": [
            {
                "edge_id": 4,
                "source": 0,
                "target": 1,
                "status": "success",
                "reoptimized": {"max_congestion": 300.0},
            },
            {
                "edge_id": 10,
                "source": 1,
                "target": 2,
                "status": "disconnected",
                "reoptimized": None,
            },
            {
                "edge_id": 14,
                "source": 2,
                "target": 3,
                "status": "success",
                "reoptimized": {"max_congestion": 240.0},
            },
        ],
    }


def test_ranking_separates_disconnected_failures():
    ranking = composite._rank_failure_outcomes(
        sample_n1_report()
    )

    assert ranking["tested_links"] == 3
    assert ranking["successful"] == 2
    assert ranking["disconnected"] == 1

    assert [
               item["edge_id"]
               for item in ranking["worst_successful_failures"]
           ] == [4, 14]

    disconnected = ranking["disconnected_failures"][0]

    assert disconnected["edge_id"] == 10
    assert disconnected["post_failure_max_congestion"] is None
    assert disconnected["relative_congestion_change"] is None

    assert ranking["worst_successful_failures"][0][
               "relative_congestion_change"
           ] == pytest.approx(0.5)


def test_composite_uses_exact_effective_workspace(
        monkeypatch,
        tmp_path,
):
    calls = []

    async def fake_dispatch(**kwargs):
        calls.append(kwargs)

        kind = kwargs["plan"].kind

        if kind == "traffic_growth":
            return {
                "network_id": kwargs["network_id"],
                "analysis_id": "traffic-001",
                "analysis_type": kind,
                "status": "completed",
                "summary": {},
                "service_result": {
                    "summary": {
                        "modified_demand_volume": 1560.0,
                    }
                },
            }

        return {
            "network_id": kwargs["network_id"],
            "analysis_id": "n1-001",
            "analysis_type": kind,
            "status": "completed",
            "summary": {
                "tested_links": 3,
                "successful": 2,
                "disconnected": 1,
                "execution_failed": 0,
            },
        }

    monkeypatch.setattr(
        composite,
        "run_analysis_plan",
        fake_dispatch,
    )

    def fake_loader(**kwargs):
        assert kwargs["network_id"] == "test-network"
        assert kwargs["analysis_id"] == "n1-001"
        assert kwargs["effective_workspace"].demand_scale == (
            pytest.approx(1.56)
        )
        return sample_n1_report()

    monkeypatch.setattr(
        composite,
        "_load_verified_n1_report",
        fake_loader,
    )

    workspace = WorkspaceSnapshot(demand_scale=1.2)
    original = deepcopy(workspace.model_dump(mode="json"))

    result = asyncio.run(
        composite.execute_traffic_n1_investigation(
            network_id="test-network",
            plan=TrafficN1InvestigationPlan(
                traffic_multiplier=1.3
            ),
            workspace=workspace,
            reports_root=tmp_path,
        )
    )

    assert result["status"] == "completed"
    assert result["traffic_analysis_id"] == "traffic-001"
    assert result["n1_analysis_id"] == "n1-001"
    assert result["effective_demand_scale"] == pytest.approx(1.56)

    assert calls[0]["plan"].kind == "traffic_growth"
    assert calls[0]["workspace"].demand_scale == pytest.approx(1.2)

    assert calls[1]["plan"].kind == "n1_resilience"
    assert calls[1]["workspace"].demand_scale == pytest.approx(1.56)

    assert workspace.model_dump(mode="json") == original

    assert result["n1_resilient"] is False
    assert result["resilience"]["disconnected"] == 1


def test_second_analysis_failure_retains_traffic_evidence(
        monkeypatch,
        tmp_path,
):
    async def fake_dispatch(**kwargs):
        if kwargs["plan"].kind == "traffic_growth":
            return {
                "network_id": kwargs["network_id"],
                "analysis_id": "traffic-002",
                "analysis_type": "traffic_growth",
                "status": "completed",
                "summary": {},
            }

        raise RuntimeError("Simulated N-1 failure")

    monkeypatch.setattr(
        composite,
        "run_analysis_plan",
        fake_dispatch,
    )

    result = asyncio.run(
        composite.execute_traffic_n1_investigation(
            network_id="test-network",
            plan=TrafficN1InvestigationPlan(
                traffic_multiplier=1.3
            ),
            workspace=WorkspaceSnapshot(),
            reports_root=tmp_path,
        )
    )

    assert result["status"] == "partial"
    assert result["traffic_analysis_id"] == "traffic-002"
    assert result["n1_analysis_id"] is None
    assert result["resilience"] is None
    assert result["n1_resilient"] is None


def test_composite_rejects_seed_mismatch(
        monkeypatch,
        tmp_path,
):
    async def never_execute(**kwargs):
        pytest.fail("No analysis should execute")

    monkeypatch.setattr(
        composite,
        "run_analysis_plan",
        never_execute,
    )

    with pytest.raises(ValueError, match="seed 42"):
        asyncio.run(
            composite.execute_traffic_n1_investigation(
                network_id="test-network",
                plan=TrafficN1InvestigationPlan(
                    traffic_multiplier=1.3,
                    seed=17,
                ),
                workspace=WorkspaceSnapshot(),
                reports_root=tmp_path,
            )
        )


def test_mismatched_demand_volumes_invalidate_comparison(
        monkeypatch,
        tmp_path,
):
    async def fake_dispatch(**kwargs):
        kind = kwargs["plan"].kind

        return {
            "network_id": kwargs["network_id"],
            "analysis_id": (
                "traffic-003"
                if kind == "traffic_growth"
                else "n1-003"
            ),
            "analysis_type": kind,
            "status": "completed",
            "summary": {},
            "service_result": (
                {
                    "summary": {
                        "modified_demand_volume": 2000.0,
                    }
                }
                if kind == "traffic_growth"
                else {}
            ),
        }

    monkeypatch.setattr(
        composite,
        "run_analysis_plan",
        fake_dispatch,
    )
    monkeypatch.setattr(
        composite,
        "_load_verified_n1_report",
        lambda **kwargs: sample_n1_report(),
    )

    result = asyncio.run(
        composite.execute_traffic_n1_investigation(
            network_id="test-network",
            plan=TrafficN1InvestigationPlan(
                traffic_multiplier=1.3
            ),
            workspace=WorkspaceSnapshot(),
            reports_root=tmp_path,
        )
    )

    assert result["status"] == "partial"
    assert result["traffic_analysis_id"] == "traffic-003"
    assert result["n1_analysis_id"] == "n1-003"
    assert result["resilience"] is None
