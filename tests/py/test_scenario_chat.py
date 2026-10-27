import json
from types import SimpleNamespace
from unittest.mock import Mock

import pytest

from e_routing.assistant.scenario_chat import answer_scenario_question


def scenario_result():
    return {
        "scenario_id": "test-123",
        "status": "completed",
        "execution_mode": "persisted_routing_scheme",
        "scenario_type": "baseline",
        "demand_model": "uniform",
        "demand_scale": 1.0,
        "capacity_upgrades": [],
        "metrics": {
            "congestion": 1.5,
            "links": [
                {"edge_id": 1, "utilization": 1.5, "overloaded": True}
            ],
            "overloaded_links": [
                {"edge_id": 1, "utilization": 1.5, "overloaded": True}
            ],
            "overloaded_count": 1,
            "layer1": None,
        },
        "limitations": ["No live telemetry."],
    }


def model_response(output, output_text=""):
    return SimpleNamespace(
        status="completed",
        output=output,
        output_text=output_text,
    )


def tool_call(name="get_scenario_summary", arguments="{}"):
    return SimpleNamespace(
        type="function_call",
        name=name,
        arguments=arguments,
        call_id="call-123",
        model_dump=lambda **kwargs: {
            "type": "function_call",
            "name": name,
            "arguments": arguments,
            "call_id": "call-123",
        },
    )


def test_evidence_required_before_answer():
    client = Mock()

    client.responses.create.side_effect = [
        model_response([tool_call()]),
        model_response([], "Congestion is 150% [E1]."),
    ]

    result = answer_scenario_question(
        scenario_result(), "What is the congestion?", client
    )

    assert result["scenario_id"] == "test-123"
    assert result["evidence"][0]["id"] == "E1"
    assert "150%" in result["answer"]

    assert client.responses.create.call_count == 2
    assert client.responses.create.call_args_list[0].kwargs["tool_choice"] == "required"
    assert client.responses.create.call_args_list[1].kwargs["tool_choice"] == "none"


def test_rejects_answer_without_tool_evidence():
    client = Mock()
    client.responses.create.return_value = model_response(
        [], "The network is safe."
    )

    with pytest.raises(RuntimeError, match="evidence"):
        answer_scenario_question(
            scenario_result(), "Is the network safe?", client
        )


def test_rejects_unauthorized_tool():
    client = Mock()
    client.responses.create.return_value = model_response(
        [tool_call(name="execute_shell")]
    )

    with pytest.raises(ValueError, match="Unauthorized"):
        answer_scenario_question(
            scenario_result(), "Run a command", client
        )


def test_rejects_unexpected_tool_arguments():
    client = Mock()
    client.responses.create.return_value = model_response(
        [tool_call(arguments=json.dumps({"scenario_id": "other"}))]
    )

    with pytest.raises(ValueError, match="arguments"):
        answer_scenario_question(
            scenario_result(), "Inspect another scenario", client
        )


def test_rejects_invalid_evidence_before_model_call():
    client = Mock()
    result = scenario_result()
    result["metrics"]["congestion"] = float("nan")

    with pytest.raises(ValueError):
        answer_scenario_question(result, "Analyze congestion", client)

    client.responses.create.assert_not_called()