
"""Regression tests for the extracted Traffic Growth service."""

from __future__ import annotations

import asyncio
import json
from pathlib import Path
from types import SimpleNamespace

import pytest
from fastapi import HTTPException

import e_routing.api.traffic_service as service
from e_routing.api.traffic_models import TrafficRequest


def sample_report(scale: float) -> dict:
    link = {
        "edge_id": 4,
        "capacity": 30.0,
        "load": 15.0,
        "utilization": 0.5,
    }

    return {
        "scenario_kind": "traffic_growth",
        "execution_mode": "fixed_routing",
        "demand_scale": scale,
        "requested_global_scale": scale,
        "absolute_demand_change": 30.0,
        "absolute_congestion_change": 0.2,
        "relative_demand_change": 0.3,
        "relative_congestion_change": 0.2,
        "commodity_overrides": [
            {
                "source": 0,
                "target": 5,
                "factor": 2.0,
            }
        ],
        "baseline": {
            "total_demand_volume": 100.0,
            "max_congestion": 1.0,
            "average_congestion": 0.4,
            "links": [dict(link)],
        },
        "modified": {
            "total_demand_volume": 130.0,
            "max_congestion": 1.2,
            "average_congestion": 0.5,
            "links": [dict(link)],
        },
    }


def configure_mock_network(monkeypatch, tmp_path):
    executable = tmp_path / "mock-engine"
    graph = tmp_path / "mock-network.lgf"
    routing_table = tmp_path / "routing.bin"

    for path in (executable, graph, routing_table):
        path.touch()

    trusted = SimpleNamespace(
        executable=executable,
        graph=graph,
        routing_table=routing_table,
        solver_label="expander",
    )

    monkeypatch.setattr(
        service,
        "load_network_registry",
        lambda: {"test-network": object()},
    )

    monkeypatch.setattr(
        service,
        "trusted_run_for_network",
        lambda network: trusted,
    )

    return trusted


def test_traffic_execution_preserves_parameters(
        monkeypatch,
        tmp_path,
):
    trusted = configure_mock_network(monkeypatch, tmp_path)
    commands = []

    def fake_run(command, **kwargs):
        commands.append(command)

        pending_path = Path(
            command[command.index("--output") + 1]
        )

        pending_path.write_text(
            json.dumps(sample_report(1.3)),
            encoding="utf-8",
        )

        return SimpleNamespace(
            returncode=0,
            stdout="",
            stderr="",
        )

    monkeypatch.setattr(
        service.subprocess,
        "run",
        fake_run,
    )

    request = TrafficRequest(
        analysis_id="traffic-test-001",
        global_scale=1.3,
        capacity_upgrades=[
            {
                "edge_id": 4,
                "capacity_factor": 1.5,
            },
        ],
        commodity_overrides=[
            {
                "source": 0,
                "target": 5,
                "factor": 2.0,
            },
        ],
    )

    response = asyncio.run(
        service.execute_traffic_analysis(
            network_id="test-network",
            request=request,
            reports_root=tmp_path,
        )
    )

    assert response["status"] == "completed"
    assert response["summary"][
               "modified_max_congestion"
           ] == 1.2

    assert len(commands) == 1
    command = commands[0]

    assert command[:2] == [
        str(trusted.executable),
        "analyze-traffic",
    ]

    assert command[
               command.index("--demand-scale") + 1
               ] == "1.3"

    assert command[
               command.index("--load-routing-table") + 1
               ] == str(trusted.routing_table)

    assert "4:1.5" in command
    assert "0:5:2.0" in command

    folder = (
            tmp_path / "test-network" / "traffic-test-001"
    )

    assert (folder / "traffic.json").is_file()
    assert not (folder / "traffic.pending.json").exists()

    metadata = json.loads(
        (folder / "traffic_request.json").read_text()
    )

    assert metadata["global_scale"] == 1.3

    status = json.loads(
        (folder / "status.json").read_text()
    )

    assert status["status"] == "completed"

    loaded = asyncio.run(
        service.load_traffic_analysis(
            "test-network",
            "traffic-test-001",
            tmp_path,
        )
    )

    assert loaded["scenario_kind"] == "traffic_growth"


def test_duplicate_traffic_analysis_is_rejected(
        monkeypatch,
        tmp_path,
):
    configure_mock_network(monkeypatch, tmp_path)

    folder = tmp_path / "test-network" / "duplicate-001"
    folder.mkdir(parents=True)

    (folder / "traffic_request.json").write_text(
        "{}",
        encoding="utf-8",
    )

    with pytest.raises(HTTPException) as caught:
        asyncio.run(
            service.execute_traffic_analysis(
                "test-network",
                TrafficRequest(
                    analysis_id="duplicate-001"
                ),
                tmp_path,
            )
        )

    assert caught.value.status_code == 409


def test_invalid_report_scale_is_not_published(
        monkeypatch,
        tmp_path,
):
    configure_mock_network(monkeypatch, tmp_path)

    def fake_run(command, **kwargs):
        pending_path = Path(
            command[command.index("--output") + 1]
        )

        # The request uses scale 1.3, but this report says 1.2.
        pending_path.write_text(
            json.dumps(sample_report(1.2)),
            encoding="utf-8",
        )

        return SimpleNamespace(
            returncode=0,
            stdout="",
            stderr="",
        )

    monkeypatch.setattr(
        service.subprocess,
        "run",
        fake_run,
    )

    with pytest.raises(HTTPException) as caught:
        asyncio.run(
            service.execute_traffic_analysis(
                "test-network",
                TrafficRequest(
                    analysis_id="invalid-scale",
                    global_scale=1.3,
                ),
                tmp_path,
            )
        )

    assert caught.value.status_code == 502

    folder = tmp_path / "test-network" / "invalid-scale"

    assert not (folder / "traffic.json").exists()

    status = json.loads(
        (folder / "status.json").read_text()
    )
    assert status["status"] == "failed"
