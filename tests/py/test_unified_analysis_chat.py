
"""Tests for Jeylan's four-operation analysis interface."""

from __future__ import annotations

import asyncio
import json
from types import SimpleNamespace

import pytest

import e_routing.assistant.unified_analysis_chat as chat
from e_routing.api.workspace_models import WorkspaceSnapshot


class FakeResponses:
    def __init__(
            self,
            operation: str,
            arguments: dict,
            explanation_fails: bool = False,
    ):
        self.operation = operation
        self.arguments = arguments
        self.explanation_fails = explanation_fails
        self.calls = 0

    def create(self, **kwargs):
        self.calls += 1

        if self.calls == 1:
            return SimpleNamespace(
                status="completed",
                output=[
                    SimpleNamespace(
                        type="function_call",
                        name=self.operation,
                        arguments=json.dumps(self.arguments),
                    )
                ],
                output_text="",
            )

        if self.explanation_fails:
            raise ConnectionError("Simulated OpenAI failure")

        return SimpleNamespace(
            status="completed",
            output=[],
            output_text="Verified analysis completed.",
        )


def run_chat(
        *,
        question: str,
        operation: str,
        arguments: dict,
        monkeypatch,
        tmp_path,
        explanation_fails: bool = False,
):
    captured = []

    async def fake_dispatch(**kwargs):
        captured.append(kwargs)

        plan = kwargs["plan"]

        return {
            "network_id": kwargs["network_id"],
            "analysis_id": "jeylan-test-001",
            "analysis_type": plan.kind,
            "status": "completed",
            "plan": plan.model_dump(mode="json"),
            "summary": {
                "baseline_max_congestion": 10.0,
                "modified_max_congestion": 12.0,
            },
            "service_result": {},
        }

    monkeypatch.setattr(
        chat,
        "run_analysis_plan",
        fake_dispatch,
    )

    client = SimpleNamespace(
        responses=FakeResponses(
            operation,
            arguments,
            explanation_fails,
        )
    )

    workspace = WorkspaceSnapshot()
    original = workspace.model_dump(mode="json")

    result = asyncio.run(
        chat.answer_unified_analysis(
            question=question,
            network_id="test-network",
            workspace=workspace,
            client=client,
            model="test-model",
            reports_root=tmp_path,
        )
    )

    assert len(captured) == 1
    assert workspace.model_dump(mode="json") == original

    return result, captured[0]


def test_capacity_upgrade_dispatch(monkeypatch, tmp_path):
    result, call = run_chat(
        question="Upgrade edge 78 by 50%",
        operation="capacity_upgrade",
        arguments={
            "demand_model": "uniform",
            "upgrades": [
                {
                    "edge_id": 78,
                    "capacity_factor": 1.5,
                }
            ],
        },
        monkeypatch=monkeypatch,
        tmp_path=tmp_path,
    )

    assert call["plan"].kind == "capacity_upgrade"
    assert call["plan"].upgrades[0].edge_id == 78
    assert result["status"] == "completed"


def test_link_failure_dispatch(monkeypatch, tmp_path):
    result, call = run_chat(
        question="Simulate failure of link 14",
        operation="link_failure",
        arguments={
            "demand_model": "uniform",
            "failed_edge_id": 14,
        },
        monkeypatch=monkeypatch,
        tmp_path=tmp_path,
    )

    assert call["plan"].kind == "link_failure"
    assert call["plan"].failed_edge_id == 14
    assert result["analysis_type"] == "link_failure"


def test_n1_dispatch(monkeypatch, tmp_path):
    result, call = run_chat(
        question="Run N-1 resilience analysis",
        operation="n1_resilience",
        arguments={
            "demand_model": "uniform",
            "execution_mode": "fast",
            "seed": 42,
        },
        monkeypatch=monkeypatch,
        tmp_path=tmp_path,
    )

    assert call["plan"].kind == "n1_resilience"
    assert call["plan"].execution_mode == "fast"
    assert result["status"] == "completed"


def test_traffic_growth_dispatch(monkeypatch, tmp_path):
    result, call = run_chat(
        question="Simulate 30% more traffic",
        operation="traffic_growth",
        arguments={
            "demand_model": "uniform",
            "traffic_multiplier": 1.3,
        },
        monkeypatch=monkeypatch,
        tmp_path=tmp_path,
    )

    assert call["plan"].kind == "traffic_growth"
    assert call["plan"].traffic_multiplier == 1.3
    assert result["analysis_type"] == "traffic_growth"


def test_invented_edge_is_rejected(monkeypatch, tmp_path):
    async def unexpected_dispatch(**kwargs):
        pytest.fail("Execution must not occur")

    monkeypatch.setattr(
        chat,
        "run_analysis_plan",
        unexpected_dispatch,
    )

    client = SimpleNamespace(
        responses=FakeResponses(
            "link_failure",
            {
                "demand_model": "uniform",
                "failed_edge_id": 78,
            },
        )
    )

    with pytest.raises(ValueError, match="physical edge ID"):
        asyncio.run(
            chat.answer_unified_analysis(
                question="Simulate failure of link 14",
                network_id="test-network",
                workspace=WorkspaceSnapshot(),
                client=client,
                model="test-model",
                reports_root=tmp_path,
            )
        )



def test_multi_analysis_request_rejects_single_tool(
        monkeypatch,
        tmp_path,
):
    """Do not silently execute only part of a composite request."""

    async def unexpected_dispatch(**kwargs):
        pytest.fail(
            "Incomplete composite request must not execute"
        )

    monkeypatch.setattr(
        chat,
        "run_analysis_plan",
        unexpected_dispatch,
    )

    monkeypatch.setattr(
        chat,
        "execute_composite_analysis",
        unexpected_dispatch,
    )

    client = SimpleNamespace(
        responses=FakeResponses(
            "traffic_growth",
            {
                "demand_model": "uniform",
                "traffic_multiplier": 1.3,
            },
        )
    )

    with pytest.raises(
            ValueError,
            match="composite",
    ):
        asyncio.run(
            chat.answer_unified_analysis(
                question=(
                    "Handle 30% more traffic while maintaining "
                    "N-1 resilience"
                ),
                network_id="test-network",
                workspace=WorkspaceSnapshot(),
                client=client,
                model="test-model",
                reports_root=tmp_path,
            )
        )



def test_completed_analysis_survives_explanation_failure(
        monkeypatch,
        tmp_path,
):
    result, call = run_chat(
        question="Simulate 30% more traffic",
        operation="traffic_growth",
        arguments={
            "demand_model": "uniform",
            "traffic_multiplier": 1.3,
        },
        monkeypatch=monkeypatch,
        tmp_path=tmp_path,
        explanation_fails=True,
    )

    assert call["plan"].kind == "traffic_growth"
    assert result["status"] == "completed"
    assert result["analysis_id"] == "jeylan-test-001"
    assert "10.0" in result["answer"]
    assert "12.0" in result["answer"]
