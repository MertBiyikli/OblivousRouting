
import pytest

from e_routing.assistant.scenario_evidence import summarize_scenario


def test_summarize_completed_scenario():
    result = {
        "scenario_id": "scenario-123",
        "status": "completed",
        "execution_mode": "persisted_routing_scheme",
        "scenario_type": "capacity_planning",
        "demand_model": "uniform",
        "demand_scale": 1.3,
        "capacity_upgrades": [{"edge_id": 2, "capacity_factor": 1.5}],
        "metrics": {
            "congestion": 2.5,
            "links": [{"edge_id": 2, "utilization": 2.5, "overloaded": True}],
            "overloaded_links": [{"edge_id": 2, "utilization": 2.5, "overloaded": True}],
            "overloaded_count": 1,
            "layer1": None,
        },
        "limitations": ["No live telemetry."],
    }

    summary = summarize_scenario(result)

    assert summary["scenario_id"] == "scenario-123"
    assert summary["evidence_source"] == "persisted_routing_scenario"
    assert summary["congestion_ratio"] == 2.5
    assert summary["overloaded_count"] == 1
    assert summary["demand_scale"] == 1.3
    assert summary["capacity_upgrades"][0]["edge_id"] == 2
    assert "run_id" not in summary


def test_rejects_incomplete_scenario():
    with pytest.raises(ValueError, match="completed"):
        summarize_scenario({"scenario_id": "running", "status": "running"})


@pytest.mark.parametrize("congestion", [
    float("nan"),
    float("inf"),
    -1.0,
    True,
    None,
])
def test_rejects_invalid_congestion(congestion):
    result = {
        "scenario_id": "test",
        "status": "completed",
        "metrics": {
            "congestion": congestion,
            "links": [],
            "overloaded_links": [],
            "overloaded_count": 0,
        },
    }

    with pytest.raises(ValueError, match="congestion"):
        summarize_scenario(result)


def test_rejects_inconsistent_overloaded_count():
    result = {
        "scenario_id": "test",
        "status": "completed",
        "metrics": {
            "congestion": 1.5,
            "links": [
                {"edge_id": 1, "utilization": 1.5, "overloaded": True}
            ],
            "overloaded_links": [],
            "overloaded_count": 0,
        },
    }

    with pytest.raises(ValueError, match="count"):
        summarize_scenario(result)