
"""Jeylan experiment tool-selection regression tests."""

from __future__ import annotations

import json
from types import SimpleNamespace

import pytest

import e_routing.assistant.experiment_chat as chat
from e_routing.api.workspace_models import WorkspaceSnapshot


class FakeResponses:
    def __init__(self, tool_name, arguments):
        self.tool_name = tool_name
        self.arguments = arguments
        self.calls = 0

    def create(self, **kwargs):
        self.calls += 1

        if self.calls == 1:
            return SimpleNamespace(
                status="completed",
                output=[
                    SimpleNamespace(
                        type="function_call",
                        name=self.tool_name,
                        arguments=json.dumps(self.arguments),
                    )
                ],
                output_text="",
            )

        return SimpleNamespace(
            status="completed",
            output=[],
            output_text="Verified experiment completed.",
        )


def make_client(tool_name, arguments):
    responses = FakeResponses(tool_name, arguments)
    return SimpleNamespace(responses=responses)


def test_bottleneck_tool_uses_source_evidence(monkeypatch):
    calls = []

    def fake_investigation(**kwargs):
        calls.append(kwargs)
        return {
            "source_evidence": {
                "source_analysis_id": "source-001",
                "candidate": {
                    "edge_id": 78,
                    "utilization": 9459.748,
                },
            },
            "experiment": {
                "analysis_id": "jeylan-test-001",
                "edge_id": 78,
                "capacity_factor": 1.5,
                "baseline_max_congestion": 9459.748,
                "fixed_routing_max_congestion": 6306.499,
                "reoptimized_max_congestion": 6678.211,
            },
        }

    monkeypatch.setattr(
        chat,
        "run_bottleneck_investigation",
        fake_investigation,
    )

    client = make_client(
        "investigate_worst_bottleneck",
        {"capacity_factor": 1.5},
    )

    result = chat.answer_capacity_experiment(
        question="Find the worst bottleneck and test a 50% upgrade",
        network_id="test-network",
        workspace=WorkspaceSnapshot(),
        source_analysis_id="source-001",
        client=client,
        model="test-model",
    )

    assert len(calls) == 1
    assert calls[0]["source_analysis_id"] == "source-001"
    assert calls[0]["capacity_factor"] == 1.5

    assert result["status"] == "completed"
    assert result["experiment"]["edge_id"] == 78
    assert len(result["evidence"]) == 2


def test_bottleneck_requires_source_analysis(monkeypatch):
    def unexpected_execution(**kwargs):
        pytest.fail("Experiment must not execute")

    monkeypatch.setattr(
        chat,
        "run_bottleneck_investigation",
        unexpected_execution,
    )

    client = make_client(
        "investigate_worst_bottleneck",
        {"capacity_factor": 1.5},
    )

    with pytest.raises(
            ValueError,
            match="Run a Capacity Planning analysis first",
    ):
        chat.answer_capacity_experiment(
            question="Investigate the worst bottleneck",
            network_id="test-network",
            workspace=WorkspaceSnapshot(),
            source_analysis_id=None,
            client=client,
            model="test-model",
        )

    assert client.responses.calls == 1


def test_model_cannot_invent_explicit_edge(monkeypatch):
    def unexpected_execution(**kwargs):
        pytest.fail("Experiment must not execute")

    monkeypatch.setattr(
        chat,
        "run_capacity_experiment",
        unexpected_execution,
    )

    client = make_client(
        "run_capacity_experiment",
        {
            "edge_id": 78,
            "capacity_factor": 1.5,
            "demand_model": "uniform",
        },
    )

    with pytest.raises(
            ValueError,
            match="Specify exactly one physical edge ID",
    ):
        chat.answer_capacity_experiment(
            question="Test a 50% capacity upgrade",
            network_id="test-network",
            workspace=WorkspaceSnapshot(),
            client=client,
            model="test-model",
        )
