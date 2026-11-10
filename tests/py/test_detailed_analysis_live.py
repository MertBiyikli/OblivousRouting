
"""Live integration tests for the detailed-analysis FastAPI endpoint.

Runs the real C++ binary through FastAPI TestClient.
No subprocess mocking is used.
"""

from __future__ import annotations

import json
import os
from pathlib import Path
from types import SimpleNamespace

import pytest
from fastapi import FastAPI
from fastapi.testclient import TestClient

from e_routing.api import detailed_analysis_api
from e_routing.api.detailed_analysis_api import (
    create_detailed_analysis_router,
)


ROOT = Path(__file__).resolve().parents[2]

GRAPH = (
        ROOT
        / "experiments/datasets/small/Backbone/1221.lgf"
)

DEFAULT_ENGINE = (
        ROOT
        / "build/mac-debug-tests-no-openmp/oblivious_routing"
)

NETWORK_ID = "backbone-1221"


@pytest.fixture
def live_client(tmp_path, monkeypatch):
    engine = Path(
        os.environ.get("ER_ENGINE", str(DEFAULT_ENGINE))
    ).resolve()

    if not engine.is_file():
        pytest.skip(
            f"C++ engine not found: {engine}"
        )

    if not GRAPH.is_file():
        pytest.skip(
            f"Test graph not found: {GRAPH}"
        )

    monkeypatch.setenv("ER_ENGINE", str(engine))

    # Keep the integration test independent of the user's
    # persisted network registry, but execute the real graph.
    network = SimpleNamespace(
        solver_label="electrical",
        graph=GRAPH,
        verify=lambda: None,
    )

    monkeypatch.setattr(
        detailed_analysis_api,
        "load_network_registry",
        lambda: {NETWORK_ID: network},
    )

    reports = tmp_path / "results" / "scenarios"

    app = FastAPI()
    app.include_router(
        create_detailed_analysis_router(
            reports_root=reports
        )
    )

    return TestClient(app), reports


def post_analysis(client, request):
    return client.post(
        f"/api/v1/networks/{NETWORK_ID}"
        "/scenarios/detailed-analysis",
        json=request,
    )


def report_for(reports, scenario_id):
    path = (
            reports
            / "detailed"
            / NETWORK_ID
            / scenario_id
            / "report.json"
    )

    assert path.is_file(), (
        f"C++ report was not created: {path}"
    )

    report = json.loads(
        path.read_text(encoding="utf-8")
    )

    assert isinstance(report, dict)
    return report


def test_live_capacity_comparison(live_client):
    client, reports = live_client

    response = post_analysis(
        client,
        {
            "scenario_id": "live-capacity-001",
            "scenario_kind": "capacity_upgrade",
            "demand_model": "uniform",
            "workspace": {
                "demand_scale": 1.0,
                "capacity_upgrades": [],
                "commodity_overrides": [],
            },
            "capacity_upgrades": [
                {
                    "edge_id": 4,
                    "capacity_factor": 1.5,
                }
            ],
            "link_failures": [],
        },
    )

    assert response.status_code == 200, response.text

    payload = response.json()

    assert payload["status"] == "completed"
    assert payload["report_available"] is True
    assert "analysis" in payload

    report = report_for(
        reports,
        "live-capacity-001",
    )

    assert report.get("status", "completed") == "completed"

    # The saved request must preserve workspace provenance.
    request_path = (
            reports
            / "detailed"
            / NETWORK_ID
            / "live-capacity-001"
            / "request.json"
    )

    request_data = json.loads(
        request_path.read_text(encoding="utf-8")
    )

    assert request_data["workspace"]["demand_scale"] == 1.0
    assert request_data["capacity_upgrades"][0]["edge_id"] == 4


def test_live_connected_link_failure(live_client):
    client, reports = live_client

    response = post_analysis(
        client,
        {
            "scenario_id": "live-failure-001",
            "scenario_kind": "link_failure",
            "demand_model": "uniform",
            "workspace": {
                "demand_scale": 1.0,
                "capacity_upgrades": [],
                "commodity_overrides": [],
            },
            "link_failures": [14],
        },
    )

    assert response.status_code == 200, response.text

    payload = response.json()

    assert payload["status"] == "completed"
    assert payload["report_available"] is True

    report = report_for(
        reports,
        "live-failure-001",
    )

    assert report.get("status", "completed") == "completed"


def test_live_modified_workspace_failure(live_client):
    client, reports = live_client

    response = post_analysis(
        client,
        {
            "scenario_id": "live-workspace-001",
            "scenario_kind": "link_failure",
            "demand_model": "uniform",
            "workspace": {
                "demand_scale": 1.3,
                "capacity_upgrades": [
                    {
                        "edge_id": 4,
                        "factor": 1.5,
                    }
                ],
                "commodity_overrides": [],
            },
            "link_failures": [14],
        },
    )

    assert response.status_code == 200, response.text

    payload = response.json()

    assert payload["status"] == "completed"

    report = report_for(
        reports,
        "live-workspace-001",
    )

    assert report.get("status", "completed") == "completed"

    request_path = (
            reports
            / "detailed"
            / NETWORK_ID
            / "live-workspace-001"
            / "request.json"
    )

    request_data = json.loads(
        request_path.read_text(encoding="utf-8")
    )

    assert request_data["workspace"]["demand_scale"] == 1.3

    assert (
            request_data["workspace"]["capacity_upgrades"][0]["edge_id"]
            == 4
    )


def test_live_disconnected_link_failure(live_client):
    edge = os.environ.get("ER_TEST_DISCONNECT_EDGE")

    if edge is None:
        pytest.skip(
            "Set ER_TEST_DISCONNECT_EDGE to a known "
            "disconnecting physical-link ID"
        )

    client, reports = live_client

    response = post_analysis(
        client,
        {
            "scenario_id": "live-disconnected-001",
            "scenario_kind": "link_failure",
            "demand_model": "uniform",
            "workspace": {
                "demand_scale": 1.0,
                "capacity_upgrades": [],
                "commodity_overrides": [],
            },
            "link_failures": [int(edge)],
        },
    )

    assert response.status_code == 200, response.text

    payload = response.json()

    assert payload["status"] == "disconnected"
    assert payload["report_available"] is True
    assert payload["failed_edge_id"] == int(edge)


    report = report_for(
        reports,
        "live-disconnected-001",
    )

    assert report["schema_version"] == "1.3"
    assert report["scenario_kind"] == "link_failure"
    assert report["link_failures"] == [int(edge)]

    solvers = report["solver_results"]
    assert len(solvers) > 0

    for solver in solvers:
        assert solver["baseline"] is not None
        assert solver["reoptimized"] is None
        assert solver["fixed_routing"] is None
        assert solver["comparison"] is None

