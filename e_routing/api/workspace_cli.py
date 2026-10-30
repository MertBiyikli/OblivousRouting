
"""Deterministic serialization of workspace CLI arguments."""

from e_routing.api.workspace_models import WorkspaceSnapshot


def workspace_cli_args(
        workspace: WorkspaceSnapshot,
        *,
        include_capacity: bool = True,
        include_commodities: bool = True,
) -> list[str]:
    workspace = workspace.normalized()

    args = [
        "--demand-scale",
        str(workspace.demand_scale),
    ]

    if include_capacity:
        for upgrade in workspace.capacity_upgrades:
            args.extend([
                "--capacity-upgrade",
                f"{upgrade.edge_id}:{upgrade.factor}",
            ])
    elif workspace.capacity_upgrades:
        raise ValueError(
            "This execution path cannot apply capacity upgrades"
        )

    if include_commodities:
        if workspace.commodity_overrides:
            args.append("--commodity-scale")

            args.extend(
                f"{item.source}:{item.target}:{item.factor}"
                for item in workspace.commodity_overrides
            )
    elif workspace.commodity_overrides:
        raise ValueError(
            "This execution path cannot apply commodity overrides"
        )

    return args
