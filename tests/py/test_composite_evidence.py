
import copy

import pytest

from e_routing.api.workspace_models import WorkspaceSnapshot
from e_routing.assistant.composite_evidence import (
    compare_n1_capacity,
)


def make_report(values, statuses=None):
    statuses = statuses or ["success"] * len(values)
    outcomes = []

    for edge_id, (value, status) in enumerate(
            zip(values, statuses)
    ):
        outcomes.append({
            "edge_id": edge_id,
            "source": edge_id,
            "target": edge_id + 1,
            "status": status,
            "reoptimized": (
                {"max_congestion": value}
                if status == "success"
                else None
            ),
        })

    return {
        "baseline": {
            "max_congestion": 100.0,
            "total_demand_volume": 2000.0,
        },
        "summary": {"tested_links": len(outcomes)},
        "outcomes": outcomes,
    }


def compare(before, after):
    return compare_n1_capacity(
        before_report=before,
        after_report=after,
        before_workspace=WorkspaceSnapshot(
            demand_scale=1.3,
        ),
        after_workspace=WorkspaceSnapshot(
            demand_scale=1.3,
            capacity_upgrades=[
                {"edge_id": 4, "factor": 1.5}
            ],
        ),
        before_analysis_id="before",
        after_analysis_id="after",
    )


def test_improvement_and_regression():
    result = compare(
        make_report([100.0, 100.0]),
        make_report([80.0, 120.0]),
    )

    assert result["counts"]["improved"] == 1
    assert result["counts"]["regressed"] == 1


def test_disconnection_recovery():
    result = compare(
        make_report(
            [None],
            ["disconnected"],
        ),
        make_report([90.0]),
    )

    assert result["counts"]["became_connected"] == 1


def test_execution_failure_not_connectivity_improvement():
    result = compare(
        make_report(
            [None],
            ["execution_failed"],
        ),
        make_report([90.0]),
    )

    assert result["counts"]["execution_status_changed"] == 1
    assert result["counts"]["became_connected"] == 0


def test_equal_congestion_is_unchanged():
    result = compare(
        make_report([100.0]),
        make_report([100.0 + 1e-11]),
    )

    assert result["counts"]["unchanged"] == 1


def test_different_edge_sets_rejected():
    with pytest.raises(ValueError, match="edge sets"):
        compare(
            make_report([100.0]),
            make_report([90.0, 95.0]),
        )


def test_different_demand_volume_rejected():
    after = make_report([90.0])
    after["baseline"]["total_demand_volume"] = 2500.0

    with pytest.raises(ValueError, match="Demand volumes"):
        compare(make_report([100.0]), after)


def test_identical_capacities_rejected():
    workspace = WorkspaceSnapshot(demand_scale=1.3)

    with pytest.raises(ValueError, match="identical"):
        compare_n1_capacity(
            before_report=make_report([100.0]),
            after_report=make_report([90.0]),
            before_workspace=workspace,
            after_workspace=workspace,
            before_analysis_id="before",
            after_analysis_id="after",
        )
