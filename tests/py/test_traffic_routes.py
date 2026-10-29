
import json
import subprocess
from pathlib import Path
from types import SimpleNamespace

import pytest
from fastapi import FastAPI
from fastapi.testclient import TestClient

from e_routing.api import traffic_routes


NETWORK = "backbone-1221"
BASE = f"/api/v1/networks/{NETWORK}/analyses/traffic"


def sample_report(scale=1.2):
    baseline = {
        "max_congestion": 2.0,
        "average_congestion": 2.0,
        "total_demand_volume": 100.0,
        "links": [{
            "edge_id": 0,
            "source": 0,
            "target": 1,
            "capacity": 50.0,
            "load": 100.0,
            "utilization": 2.0,
        }],
    }

    modified = {
        "max_congestion": 2.0 * scale,
        "average_congestion": 2.0 * scale,
        "total_demand_volume": 100.0 * scale,
        "links": [{
            "edge_id": 0,
            "source": 0,
            "target": 1,
            "capacity": 50.0,
            "load": 100.0 * scale,
            "utilization": 2.0 * scale,
        }],
    }

    return {
        "schema_version": "1.0",
        "scenario_kind": "traffic_growth",
        "execution_mode": "fixed_routing",
        "graph_fingerprint": "test-fingerprint",
        "baseline": baseline,
        "modified": modified,
        "demand_scale": scale,
        "requested_global_scale": scale,
        "commodity_overrides": [],
        "absolute_demand_change": 100.0 * (scale - 1),
        "relative_demand_change": scale - 1,
        "absolute_congestion_change": 2.0 * (scale - 1),
        "relative_congestion_change": scale - 1,
    }


@pytest.fixture
def api(tmp_path, monkeypatch):
    # No real registered network or executable is used.
    trusted = SimpleNamespace(
        executable=Path("/fake/engine"),
        graph=Path("/fake/graph.lgf"),
        routing_table=Path("/fake/routing.bin"),
        solver_label="electrical",
    )

    monkeypatch.setattr(
        traffic_routes,
        "load_network_registry",
        lambda: {NETWORK: object()},
    )

    monkeypatch.setattr(
        traffic_routes,
        "trusted_run_for_network",
        lambda network: trusted,
    )

    root = tmp_path / "reports"

    app = FastAPI()
    app.include_router(
        traffic_routes.create_traffic_router(
            reports_root=root,
            timeout_seconds=1,
        )
    )

    with TestClient(app) as client:
        yield client, root


def payload(analysis_id="test-001", scale=1.2):
    return {
        "analysis_id": analysis_id,
        "demand_model": "uniform",
        "seed": 42,
        "global_scale": scale,
        "commodity_overrides": [],
    }


def folder(root, analysis_id):
    return root / NETWORK / analysis_id


def status(root, analysis_id):
    return json.loads(
        (folder(root, analysis_id) / "status.json").read_text()
    )


def get_report(client, analysis_id):
    return client.get(f"{BASE}/{analysis_id}")


def mock_execution(monkeypatch, behavior):
    def run(command, **kwargs):
        assert command[1] == "analyze-traffic"
        assert "--load-routing-table" in command
        assert kwargs["timeout"] == 1

        output = Path(command[command.index("--output") + 1])
        scale = float(
            command[command.index("--demand-scale") + 1]
        )

        if behavior == "timeout":
            raise subprocess.TimeoutExpired(command, 1)

        if behavior == "engine_error":
            return subprocess.CompletedProcess(
                command, 1, stdout="", stderr="Simulated failure"
            )

        output.parent.mkdir(parents=True, exist_ok=True)

        if behavior == "invalid_json":
            output.write_text("{invalid-json")
        else:
            output.write_text(
                json.dumps(sample_report(scale))
            )

        return subprocess.CompletedProcess(
            command, 0, stdout="OK", stderr=""
        )

    monkeypatch.setattr(traffic_routes.subprocess, "run", run)


def test_successful_analysis(api, monkeypatch):
    client, root = api
    mock_execution(monkeypatch, "success")

    response = client.post(BASE, json=payload())

    assert response.status_code == 201, response.text
    assert response.json()["status"] == "completed"
    assert status(root, "test-001")["status"] == "completed"

    result = get_report(client, "test-001")

    assert result.status_code == 200
    assert result.json()["requested_global_scale"] == 1.2
    assert result.json()["modified"]["max_congestion"] == 2.4


@pytest.mark.parametrize(
    "behavior,expected_http",
    [
        ("engine_error", 502),
        ("timeout", 504),
        ("invalid_json", 502),
    ],
)
def test_failed_execution(
        api, monkeypatch, behavior, expected_http
):
    client, root = api
    analysis_id = f"failed-{behavior}"

    mock_execution(monkeypatch, behavior)

    response = client.post(
        BASE, json=payload(analysis_id)
    )

    assert response.status_code == expected_http
    assert status(root, analysis_id)["status"] == "failed"

    # The original request remains available for diagnosis.
    request_file = folder(
        root, analysis_id
    ) / "traffic_request.json"

    assert request_file.is_file()
    assert json.loads(
        request_file.read_text()
    )["analysis_id"] == analysis_id

    # No completed report may be published.
    assert not (
            folder(root, analysis_id) / "traffic.json"
    ).exists()

    result = get_report(client, analysis_id)

    assert result.status_code == 200
    assert result.json()["status"] == "failed"


def test_duplicate_does_not_overwrite(api, monkeypatch):
    client, root = api
    mock_execution(monkeypatch, "success")

    first = client.post(BASE, json=payload(scale=1.2))
    assert first.status_code == 201

    duplicate = client.post(
        BASE, json=payload(scale=1.5)
    )

    assert duplicate.status_code == 409

    report = get_report(client, "test-001")
    assert report.status_code == 200
    assert report.json()["requested_global_scale"] == 1.2
    assert status(root, "test-001")["status"] == "completed"


def test_invalid_input_creates_no_reservation(api):
    client, root = api

    bad = payload("invalid-input")
    bad["global_scale"] = -2.0

    response = client.post(BASE, json=bad)

    assert response.status_code == 422
    assert not folder(root, "invalid-input").exists()
