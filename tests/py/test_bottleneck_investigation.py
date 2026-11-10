
"""Tests for deterministic evidence-driven bottleneck selection."""

from __future__ import annotations

from types import SimpleNamespace

import pytest

import e_routing.assistant.bottleneck_investigation as investigation
from e_routing.api.workspace_models import WorkspaceSnapshot


def sample_report():
    return {
        "scenario_kind": "capacity_upgrade",
        "solver_results": [{
            "baseline": {
                "links": [
                    {
                        "edge_id": 4,
                        "source": 0,
                        "target": 1,
                        "load": 40.0,
                        "capacity": 10.0,
                        "utilization": 4.0,
                    },
                    {
                        "edge_id": 78,
                        "source": 3,
                        "target": 9,
                        "load": 100.0,
                        "capacity": 10.0,
                        "utilization": 10.0,
                    },
                    {
                        "edge_id": 12,
                        "source": 1,
                        "target": 2,
                        "load": 60.0,
                        "capacity": 10.0,
                        "utilization": 6.0,
                    },
                ],
            },
        }],
    }


def test_selects_worst_baseline_bottleneck():
    selected = investigation.select_worst_baseline_link(
        sample_report()
    )

    assert selected["edge_id"] == 78
    assert selected["utilization"] == 10.0
    assert selected["capacity"] == 10.0


def test_tie_breaks_by_smallest_edge_id():
    report = sample_report()

    report["solver_results"][0]["baseline"]["links"][0][
        "utilization"
    ] = 10.0

    selected = investigation.select_worst_baseline_link(
        report
    )

    assert selected["edge_id"] == 4


@pytest.mark.parametrize(
    "invalid_value",
    [float("nan"), float("inf"), -1.0],
)
def test_rejects_invalid_utilization(invalid_value):
    report = sample_report()

    report["solver_results"][0]["baseline"]["links"][0][
        "utilization"
    ] = invalid_value

    with pytest.raises(ValueError):
        investigation.select_worst_baseline_link(report)


def test_investigation_runs_one_selected_experiment(
        monkeypatch,
        tmp_path,
):
    calls = []

    evidence = {
        "source_analysis_id": "source-001",
        "network_id": "test-network",
        "demand_model": "uniform",
        "ranking": "highest_baseline_utilization",
        "candidate": {
            "edge_id": 78,
            "utilization": 10.0,
        },
        "evaluated_links": 3,
    }

    monkeypatch.setattr(
        investigation,
        "load_bottleneck_evidence",
        lambda **kwargs: evidence,
    )

    def fake_execute(**kwargs):
        calls.append(kwargs)
        return SimpleNamespace(
            model_dump=lambda **options: {
                "analysis_id": "jeylan-test-001",
                "edge_id": 78,
                "capacity_factor": 1.5,
                "baseline_max_congestion": 10.0,
                "fixed_routing_max_congestion": 7.0,
                "reoptimized_max_congestion": 8.0,
            }
        )

    monkeypatch.setattr(
        investigation,
        "run_capacity_experiment",
        fake_execute,
    )

    workspace = WorkspaceSnapshot()
    original = workspace.model_dump(mode="json")

    response = investigation.run_bottleneck_investigation(
        network_id="test-network",
        source_analysis_id="source-001",
        workspace=workspace,
        capacity_factor=1.5,
        reports_root=tmp_path,
    )

    assert len(calls) == 1
    assert calls[0]["plan"].edge_id == 78
    assert calls[0]["plan"].capacity_factor == 1.5
    assert calls[0]["plan"].demand_model == "uniform"

    assert response["experiment"]["edge_id"] == 78
    assert response["source_evidence"]["candidate"][
               "utilization"
           ] == 10.0

    assert workspace.model_dump(mode="json") == original


def test_existing_workspace_upgrades_are_rejected(
        monkeypatch,
):
    def must_not_load(**kwargs):
        pytest.fail("Must reject before loading evidence")

    monkeypatch.setattr(
        investigation,
        "load_bottleneck_evidence",
        must_not_load,
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
        investigation.run_bottleneck_investigation(
            network_id="test-network",
            source_analysis_id="source-001",
            workspace=workspace,
        )
