
import json
from types import SimpleNamespace

import pytest

from e_routing.api.workspace_models import WorkspaceSnapshot
from e_routing.assistant import unified_analysis_chat as chat


QUESTION = (
    "Increase traffic by 30%, evaluate N-1, "
    "upgrade edge 4 by 50%, then evaluate N-1 again."
)

ARGUMENTS = {
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
                {"edge_id": 4, "capacity_factor": 1.5}
            ],
        },
        {
            "operation": "n1_resilience",
            "execution_mode": "fast",
        },
    ],
}


def tool_call(arguments=None):
    return SimpleNamespace(
        name="composite_analysis",
        arguments=json.dumps(
            arguments if arguments is not None else ARGUMENTS
        ),
        type="function_call",
    )


def response(output=None, output_text=""):
    return SimpleNamespace(
        status="completed",
        output=output or [],
        output_text=output_text,
    )


def fake_client(*responses):
    iterator = iter(responses)

    return SimpleNamespace(
        responses=SimpleNamespace(
            create=lambda **kwargs: next(iterator)
        )
    )


def fake_result(status="completed"):
    completed = [
        {
            "step_index": i,
            "operation": operation,
            "analysis_id": f"test-{i}",
            "status": "completed",
            "plan": {"kind": operation},
            "summary": {},
            "input_workspace": {},
            "output_workspace": {},
        }
        for i, operation in enumerate([
            "traffic_growth",
            "n1_resilience",
            "capacity_upgrade",
            "n1_resilience",
        ])
    ]

    return {
        "network_id": "backbone-1221",
        "kind": "composite",
        "status": status,
        "plan": {"kind": "composite"},
        "completed_steps": completed,
        "comparisons": [],
        "comparison_status": "completed",
        "effective_workspace": {},
        "workspace_modified": False,
        "failed_step_index": None,
    }


def test_explicit_capacity_edge_accepted():
    plan = chat._validated_composite_plan(
        tool_call(), QUESTION
    )

    assert plan.kind == "composite"
    assert len(plan.steps) == 4


def test_invented_edge_rejected():
    arguments = json.loads(json.dumps(ARGUMENTS))
    arguments["steps"][2]["upgrades"][0]["edge_id"] = 99

    with pytest.raises(ValueError, match="edge IDs"):
        chat._validated_composite_plan(
            tool_call(arguments),
            QUESTION,
        )


def test_multiple_capacity_events_rejected():
    arguments = json.loads(json.dumps(ARGUMENTS))

    arguments["steps"].append({
        "operation": "capacity_upgrade",
        "upgrades": [
            {"edge_id": 10, "capacity_factor": 1.5}
        ],
    })

    with pytest.raises(
            ValueError, match="invalid composite plan"
    ):
        chat._validated_composite_plan(
            tool_call(arguments),
            QUESTION + " Also upgrade edge 10.",
            )


@pytest.mark.anyio
async def test_composite_dispatch(monkeypatch):
    calls = []

    async def execute(**kwargs):
        calls.append(kwargs)
        return fake_result()

    monkeypatch.setattr(
        chat, "execute_composite_analysis", execute
    )

    client = fake_client(
        response(output=[tool_call()]),
        response(output_text="Verified analysis explanation."),
    )

    result = await chat.answer_unified_analysis(
        question=QUESTION,
        network_id="backbone-1221",
        workspace=WorkspaceSnapshot(),
        client=client,
        model="test-model",
    )

    assert len(calls) == 1
    assert result["status"] == "completed"
    assert result["analysis_type"] == "composite"
    assert len(result["completed_steps"]) == 4
    assert len(result["evidence"]) == 4
    assert result["analysis_id"] is None


@pytest.mark.anyio
async def test_partial_execution_preserved(monkeypatch):
    async def execute(**kwargs):
        result = fake_result("partial")
        result["completed_steps"] = result["completed_steps"][:2]
        result["failed_step_index"] = 2
        result["failed_operation"] = "capacity_upgrade"
        result["error"] = "Step failed"
        return result

    monkeypatch.setattr(
        chat, "execute_composite_analysis", execute
    )

    client = fake_client(
        response(output=[tool_call()]),
        response(output_text=""),
    )

    result = await chat.answer_unified_analysis(
        question=QUESTION,
        network_id="backbone-1221",
        workspace=WorkspaceSnapshot(),
        client=client,
        model="test-model",
    )

    assert result["status"] == "partial"
    assert len(result["completed_steps"]) == 2
    assert len(result["evidence"]) == 2
    assert "did not fully complete" in result["answer"]
