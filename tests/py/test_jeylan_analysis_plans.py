
"""Tests for Jeylan's unified analysis-plan contract."""

import pytest
from pydantic import TypeAdapter, ValidationError

from e_routing.assistant.analysis_plan import (
    AnalysisOperation,
)


adapter = TypeAdapter(AnalysisOperation)


@pytest.mark.parametrize(
    "payload,expected_kind",
    [
        (
                {
                    "kind": "capacity_upgrade",
                    "upgrades": [
                        {
                            "edge_id": 78,
                            "capacity_factor": 1.5,
                        }
                    ],
                },
                "capacity_upgrade",
        ),
        (
                {
                    "kind": "link_failure",
                    "failed_edge_id": 14,
                },
                "link_failure",
        ),
        (
                {
                    "kind": "n1_resilience",
                    "execution_mode": "fast",
                },
                "n1_resilience",
        ),
        (
                {
                    "kind": "traffic_growth",
                    "traffic_multiplier": 1.3,
                },
                "traffic_growth",
        ),
    ],
)
def test_supported_operations(payload, expected_kind):
    plan = adapter.validate_python(payload)
    assert plan.kind == expected_kind


@pytest.mark.parametrize(
    "payload",
    [
        {
            "kind": "link_failure",
            "failed_edge_id": -1,
        },
        {
            "kind": "traffic_growth",
            "traffic_multiplier": float("nan"),
        },
        {
            "kind": "traffic_growth",
            "traffic_multiplier": 0,
        },
        {
            "kind": "n1_resilience",
            "execution_mode": "research",
        },
        {
            "kind": "capacity_upgrade",
            "upgrades": [],
        },
        {
            "kind": "capacity_upgrade",
            "upgrades": [
                {
                    "edge_id": 4,
                    "capacity_factor": 1.5,
                },
                {
                    "edge_id": 10,
                    "capacity_factor": 1.5,
                },
            ],
        },
        {
            "kind": "unsupported_analysis",
        },
    ],
)
def test_invalid_operations_are_rejected(payload):
    with pytest.raises(ValidationError):
        adapter.validate_python(payload)
