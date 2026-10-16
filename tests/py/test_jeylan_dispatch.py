
from unittest.mock import Mock

import pytest

from e_routing.assistant.tools import ToolDispatcher
from scripts.jeylan import JeylanAssistant


def test_dispatches_capacity_planning_scenario():
    bridge = Mock()
    bridge.evaluate.return_value = {
        "scenario_id": "test-scenario",
        "status": "completed",
    }

    assistant = JeylanAssistant(
        analysis_api=Mock(),
        client=Mock(),
        run_id="test-run",
        scenario_bridge=bridge,
    )

    args = {
        "demand_model": "uniform",
        "include_layer1": False,
        "demand_scale": 1.3,
        "capacity_upgrades": [
            {"edge_id": 4, "capacity_factor": 2.0}
        ],
    }

    result = assistant._dispatch("evaluate_persisted_scenario", args)

    bridge.evaluate.assert_called_once_with(args)
    assert result["status"] == "completed"




def test_scenario_dispatch():
    bridge = Mock()
    bridge.evaluate.return_value = {"scenario_id": "scenario-1", "status": "completed"}

    dispatcher = ToolDispatcher(Mock(), "run-1", scenario_bridge=bridge)

    arguments = {
        "demand_model": "uniform",
        "include_layer1": False,
        "demand_scale": 1.3,
        "capacity_upgrades": [{"edge_id": 4, "capacity_factor": 2.0}],
    }

    result = dispatcher.dispatch("evaluate_persisted_scenario", arguments)

    assert result["status"] == "completed"
    bridge.evaluate.assert_called_once_with(arguments)


def test_unknown_tool_rejected():
    dispatcher = ToolDispatcher(Mock(), "run-1")

    with pytest.raises(ValueError, match="Unknown/unauthorized"):
        dispatcher.dispatch("execute_shell_command", {})


def test_model_cannot_select_run():
    api = Mock()
    dispatcher = ToolDispatcher(api, "trusted-run")

    dispatcher.dispatch("get_analysis_summary", {})

    api.get_analysis_summary.assert_called_once_with("trusted-run")


def test_invalid_edge_id_rejected():
    dispatcher = ToolDispatcher(Mock(), "run-1")

    with pytest.raises(ValueError):
        dispatcher.dispatch("get_failure_impact", {"failed_edge_id": -1})
