
import pytest
from pydantic import ValidationError

from e_routing.api.workspace_models import WorkspaceSnapshot


def test_empty_workspace():
    workspace = WorkspaceSnapshot()

    assert workspace.demand_scale == 1.0
    assert workspace.capacity_upgrades == []
    assert workspace.commodity_overrides == []


def test_combined_workspace():
    workspace = WorkspaceSnapshot(
        capacity_upgrades=[
            {"edge_id": 4, "factor": 1.5},
            {"edge_id": 8, "factor": 2.0},
        ],
        demand_scale=1.3,
        commodity_overrides=[
            {"source": 0, "target": 2, "factor": 1.5},
        ],
    )

    assert len(workspace.capacity_upgrades) == 2
    assert workspace.demand_scale == 1.3
    assert len(workspace.commodity_overrides) == 1


def test_duplicate_capacity_rejected():
    with pytest.raises(ValidationError):
        WorkspaceSnapshot(
            capacity_upgrades=[
                {"edge_id": 4, "factor": 1.5},
                {"edge_id": 4, "factor": 2.0},
            ]
        )


def test_invalid_demand_scale_rejected():
    with pytest.raises(ValidationError):
        WorkspaceSnapshot(demand_scale=float("inf"))


def test_duplicate_commodity_rejected():
    with pytest.raises(ValidationError):
        WorkspaceSnapshot(
            commodity_overrides=[
                {"source": 0, "target": 2, "factor": 1.5},
                {"source": 0, "target": 2, "factor": 2.0},
            ]
        )


def test_normalization():
    workspace = WorkspaceSnapshot(
        capacity_upgrades=[
            {"edge_id": 8, "factor": 2.0},
            {"edge_id": 4, "factor": 1.5},
        ],
    )

    normalized = workspace.normalized()

    assert [
               item.edge_id for item in normalized.capacity_upgrades
           ] == [4, 8]



from e_routing.api.n1_routes import N1Request


def test_n1_accepts_combined_workspace():
    request = N1Request(
        analysis_id="combined-workspace-001",
        demand_model="uniform",
        seed=42,
        execution_mode="fast",
        workspace={
            "capacity_upgrades": [
                {"edge_id": 4, "factor": 1.5}
            ],
            "demand_scale": 1.3,
            "commodity_overrides": [
                {"source": 0, "target": 2, "factor": 1.5}
            ],
        },
    )

    assert request.workspace.demand_scale == 1.3
    assert request.workspace.capacity_upgrades[0].edge_id == 4
    assert request.workspace.commodity_overrides[0].factor == 1.5
