
"""Regression tests for the extracted N-1 application service."""

from __future__ import annotations

import asyncio
import json
from pathlib import Path
from types import SimpleNamespace

import pytest
from fastapi import HTTPException

import e_routing.api.n1_service as service


def test_research_mode_is_rejected_before_execution(
        tmp_path,
        monkeypatch,
):
    def unexpected_registry():
        pytest.fail("Registry must not be read")

    monkeypatch.setattr(
        service,
        "load_network_registry",
        unexpected_registry,
    )

    request = service.N1Request(
        analysis_id="n1-research",
        execution_mode="research",
    )

    with pytest.raises(HTTPException) as caught:
        asyncio.run(
            service.execute_n1_analysis(
                network_id="test-network",
                request=request,
                reports_root=tmp_path,
            )
        )

    assert caught.value.status_code == 422
    assert "RESEARCH" in caught.value.detail


def test_n1_service_preserves_workspace_and_report(
        tmp_path,
        monkeypatch,
):
    executable = tmp_path / "mock-engine"
    executable.touch()

    graph = tmp_path / "mock-network.lgf"
    graph.touch()

    network = SimpleNamespace(
        solver_label="expander",
    )

    trusted = SimpleNamespace(
        executable=executable,
        graph=graph,
    )

    monkeypatch.setattr(
        service,
        "load_network_registry",
        lambda: {"test-network": network},
    )

    monkeypatch.setattr(
        service,
        "trusted_run_for_network",
        lambda selected_network: trusted,
    )

    # The production report validator is reused unchanged.
    # Mock only its output for this service-boundary test.
    monkeypatch.setattr(
        service,
        "validate_n1_report",
        lambda report: report,
    )

    captured_commands = []

    async def fake_subprocess(*command, **kwargs):
        captured_commands.append(list(command))

        output_index = command.index("--output") + 1
        pending_path = Path(command[output_index])

        pending_path.write_text(
            json.dumps({
                "solver": "expander_hierarchy",
                "summary": {
                    "tested_links": 3,
                    "successful": 2,
                    "disconnected": 1,
                    "execution_failed": 0,
                },
            }),
            encoding="utf-8",
        )

        class FakeProcess:
            returncode = 0

            async def communicate(self):
                return b"", b""

        return FakeProcess()

    monkeypatch.setattr(
        service.asyncio,
        "create_subprocess_exec",
        fake_subprocess,
    )

    request = service.N1Request(
        analysis_id="n1-test-001",
        demand_model="uniform",
        workspace={
            "demand_scale": 1.3,
            "capacity_upgrades": [
                {"edge_id": 4, "factor": 1.5},
            ],
            "commodity_overrides": [
                {
                    "source": 0,
                    "target": 5,
                    "factor": 2.0,
                },
            ],
        },
    )

    response = asyncio.run(
        service.execute_n1_analysis(
            network_id="test-network",
            request=request,
            reports_root=tmp_path,
        )
    )

    assert response["status"] == "completed"
    assert response["summary"]["tested_links"] == 3

    assert len(captured_commands) == 1
    command = captured_commands[0]

    assert command[0] == str(executable)
    assert command[1] == "analyze-n1"

    assert command[
               command.index("--demand-scale") + 1
               ] == "1.3"

    assert "--capacity-upgrade" in command
    assert "4:1.5" in command

    assert "--commodity-scale" in command
    assert "0:5:2.0" in command

    folder = tmp_path / "test-network" / "n1-test-001"

    assert (folder / "request.json").is_file()
    assert (folder / "n1.json").is_file()
    assert not (folder / "n1.pending.json").exists()

    metadata = json.loads(
        (folder / "request.json").read_text()
    )

    assert metadata["network_id"] == "test-network"
    assert metadata["request"]["workspace"][
               "demand_scale"
           ] == 1.3


def test_duplicate_analysis_id_is_rejected(
        tmp_path,
        monkeypatch,
):
    executable = tmp_path / "engine"
    executable.touch()

    network = SimpleNamespace(
        solver_label="electrical",
    )

    monkeypatch.setattr(
        service,
        "load_network_registry",
        lambda: {"test-network": network},
    )

    monkeypatch.setattr(
        service,
        "trusted_run_for_network",
        lambda selected_network: SimpleNamespace(
            executable=executable,
            graph=tmp_path / "graph.lgf",
        ),
    )

    folder = tmp_path / "test-network" / "duplicate-001"
    folder.mkdir(parents=True)

    (folder / "request.json").write_text(
        '{"status":"running"}',
        encoding="utf-8",
    )

    request = service.N1Request(
        analysis_id="duplicate-001",
    )

    with pytest.raises(HTTPException) as caught:
        asyncio.run(
            service.execute_n1_analysis(
                network_id="test-network",
                request=request,
                reports_root=tmp_path,
            )
        )

    assert caught.value.status_code == 409
