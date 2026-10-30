
import pytest

from e_routing.api.workspace_cli import workspace_cli_args
from e_routing.api.workspace_models import WorkspaceSnapshot


def test_empty_workspace():
    workspace = WorkspaceSnapshot()

    assert workspace_cli_args(workspace) == [
        "--demand-scale", "1.0",
    ]


def test_combined_workspace_arguments():
    workspace = WorkspaceSnapshot.model_validate({
        "demand_scale": 1.3,
        "capacity_upgrades": [
            {"edge_id": 4, "factor": 1.5},
        ],
        "commodity_overrides": [
            {
                "source": 0,
                "target": 2,
                "factor": 1.5,
            },
        ],
    })

    assert workspace_cli_args(workspace) == [
        "--demand-scale", "1.3",
        "--capacity-upgrade", "4:1.5",
        "--commodity-scale", "0:2:1.5",
    ]


def test_unsupported_capacity_is_rejected():
    workspace = WorkspaceSnapshot.model_validate({
        "capacity_upgrades": [
            {"edge_id": 4, "factor": 1.5}
        ]
    })

    with pytest.raises(ValueError, match="capacity"):
        workspace_cli_args(
            workspace,
            include_capacity=False,
        )


def test_unsupported_commodity_is_rejected():
    workspace = WorkspaceSnapshot.model_validate({
        "commodity_overrides": [
            {
                "source": 0,
                "target": 2,
                "factor": 1.5,
            }
        ]
    })

    with pytest.raises(ValueError, match="commodity"):
        workspace_cli_args(
            workspace,
            include_commodities=False,
        )
