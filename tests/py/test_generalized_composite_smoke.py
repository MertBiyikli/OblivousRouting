
"""Live generalized composite analysis smoke test."""

import asyncio
import os
from pathlib import Path


from e_routing.api.workspace_models import WorkspaceSnapshot
from e_routing.assistant.composite_plans import (
    CompositeAnalysisPlan,
)
from e_routing.assistant.composite_executor import (
    execute_composite_analysis,
)


async def main():
    network_id = os.environ.get(
        "ER_SMOKE_NETWORK_ID",
        "backbone-1221",
    )

    workspace = WorkspaceSnapshot()

    original = workspace.model_dump(mode="json")

    plan = CompositeAnalysisPlan.model_validate({
        "kind": "composite",
        "demand_model": "uniform",
        "seed": 42,
        "steps": [
            {
                "operation": "traffic_growth",
                "traffic_multiplier": 1.3,
            },
            {
                "operation": "n1_resilience",
                "execution_mode": "fast",
            },
            {
                "operation": "capacity_upgrade",
                "upgrades": [
                    {
                        "edge_id": 4,
                        "capacity_factor": 1.5,
                    }
                ],
            },
            {
                "operation": "n1_resilience",
                "execution_mode": "fast",
            },
        ],
    })

    result = await execute_composite_analysis(
        network_id=network_id,
        plan=plan,
        workspace=workspace,
        reports_root=Path("results/scenarios"),
    )

    print("Status:", result["status"])
    print("Completed steps:", len(result["completed_steps"]))
    print("Comparisons:", len(result.get("comparisons", [])))

    for step in result["completed_steps"]:
        print(
            step["operation"],
            step["analysis_id"],
        )

    assert result["status"] == "completed", result
    assert len(result["completed_steps"]) == 4
    assert result["comparison_status"] == "completed"
    assert len(result["comparisons"]) == 1

    comparison = result["comparisons"][0]
    assert comparison["comparison_type"] == "capacity_n1"

    steps = result["completed_steps"]

    assert steps[1]["input_workspace"]["demand_scale"] == 1.3
    assert steps[3]["input_workspace"]["demand_scale"] == 1.3

    assert steps[1]["input_workspace"]["capacity_upgrades"] == []
    assert steps[3]["input_workspace"]["capacity_upgrades"] == [
        {"edge_id": 4, "factor": 1.5}
    ]

    assert result["workspace_modified"] is False
    assert workspace.model_dump(mode="json") == original

    print("GENERALIZED COMPOSITE SMOKE: PASS")


if __name__ == "__main__":
    asyncio.run(main())
