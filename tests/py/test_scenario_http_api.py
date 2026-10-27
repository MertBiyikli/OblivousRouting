
"""Tests for the E-Routing HTTP API."""

import os
from pathlib import Path
from unittest.mock import patch

import pytest
from fastapi.testclient import TestClient

from e_routing.api.main import app
from e_routing.scenarios.repository import ScenarioRepository
from e_routing.scenarios.runtime import ScenarioRequest, TrustedRun

client = TestClient(app)


@pytest.fixture
def scenario_environment(tmp_path, monkeypatch):
    """Provide isolated scenario storage and trusted configuration."""
    monkeypatch.setenv("ER_RESPONSE_DIR", str(tmp_path / "responses"))

    trusted = TrustedRun(
        executable=tmp_path / "engine",
        graph=tmp_path / "graph.lgf",
        routing_table=tmp_path / "routing.bin",
        solver_label="electrical",
        workspace=tmp_path / "executions",
    )

    with patch("e_routing.api.main.load_trusted_config", return_value=trusted):
        yield trusted



def valid_payload():
    return {
        "scenario_id": "test-capacity-01",
        "scenario_type": "capacity_planning",
        "demand_model": "uniform",
        "include_layer1": False,
        "demand_scale": 1.3,
        "capacity_upgrades": [
            {"edge_id": 2, "capacity_factor": 1.5}
        ],
    }


def test_health():
    response = client.get("/health")
    assert response.status_code == 200
    assert response.json()["status"] == "ok"


def test_valid_scenario(scenario_environment):
    expected = {
        "scenario_id": "test-capacity-01",
        "status": "completed",
        "metrics": {"congestion": 12.5},
    }

    with patch("e_routing.api.main.execute", return_value=expected) as execute:
        response = client.post("/api/v1/scenarios", json=valid_payload())

    assert response.status_code == 200
    assert response.json() == expected

    request, trusted = execute.call_args.args
    assert isinstance(request, ScenarioRequest)
    assert trusted.workspace != scenario_environment.workspace
    assert trusted.workspace.parent == scenario_environment.workspace

    repository = ScenarioRepository(
        Path(os.environ["ER_RESPONSE_DIR"])
    )
    assert repository.get("test-capacity-01") == expected


def test_rejects_invalid_demand_scale():
    payload = valid_payload()
    payload["demand_scale"] = -1

    response = client.post("/api/v1/scenarios", json=payload)

    assert response.status_code == 422


def test_rejects_unknown_fields():
    payload = valid_payload()
    payload["executable"] = "/tmp/malicious"

    response = client.post("/api/v1/scenarios", json=payload)

    assert response.status_code == 422


def test_rejects_duplicate_capacity_edges():
    payload = valid_payload()
    payload["capacity_upgrades"].append(
        {"edge_id": 2, "capacity_factor": 2.0}
    )

    response = client.post("/api/v1/scenarios", json=payload)

    assert response.status_code == 422



def test_missing_server_configuration():
    with patch(
            "e_routing.api.main.load_trusted_config",
            side_effect=RuntimeError("Missing configuration")
    ):
        response = client.post("/api/v1/scenarios", json=valid_payload())

    assert response.status_code == 503
    assert response.json()["detail"] == "Scenario service unavailable"


def test_engine_execution_failure(scenario_environment):
    with patch(
            "e_routing.api.main.execute",
            side_effect=RuntimeError("Engine failed")
    ):
        response = client.post("/api/v1/scenarios", json=valid_payload())

    assert response.status_code == 502
    assert response.json()["detail"] == "Scenario execution failed"

    repository = ScenarioRepository(Path(os.environ["ER_RESPONSE_DIR"]))
    assert repository.get("test-capacity-01") is None


def test_engine_error_does_not_leak_details(scenario_environment):
    secret = "/private/internal/engine/path"

    with patch(
            "e_routing.api.main.execute",
            side_effect=RuntimeError(secret)
    ):
        response = client.post("/api/v1/scenarios", json=valid_payload())

    assert response.status_code == 502
    assert secret not in response.text



def test_duplicate_scenario_returns_conflict(scenario_environment):
    expected = {
        "scenario_id": "test-capacity-01",
        "status": "completed",
        "metrics": {"congestion": 12.5},
    }

    with patch("e_routing.api.main.execute", return_value=expected) as execute:
        first = client.post("/api/v1/scenarios", json=valid_payload())
        second = client.post("/api/v1/scenarios", json=valid_payload())

    assert first.status_code == 200
    assert second.status_code == 409
    assert second.json()["detail"] == "Scenario ID already exists"
    execute.assert_called_once()



def test_get_completed_scenario(scenario_environment):
    repository = ScenarioRepository(Path(os.environ["ER_RESPONSE_DIR"]))
    expected = {"scenario_id": "existing-01", "status": "completed"}
    repository.save(expected)

    response = client.get("/api/v1/scenarios/existing-01")

    assert response.status_code == 200
    assert response.json() == expected


def test_get_running_scenario(scenario_environment):
    repository = ScenarioRepository(Path(os.environ["ER_RESPONSE_DIR"]))
    repository.reserve("running-01")

    response = client.get("/api/v1/scenarios/running-01")

    assert response.status_code == 409
    assert response.json()["detail"] == "Scenario is still running"


@pytest.mark.parametrize("scenario_id", ["übung", "scénario", "测试"])
def test_rejects_non_ascii_scenario_id(scenario_id):
    payload = valid_payload()
    payload["scenario_id"] = scenario_id

    response = client.post("/api/v1/scenarios", json=payload)

    assert response.status_code == 422