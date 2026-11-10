
import pytest
from pydantic import ValidationError

from e_routing.api.workspace_models import WorkspaceSnapshot
from e_routing.assistant.composite_plans import (
    TrafficN1InvestigationPlan,
    prepare_traffic_n1_investigation,
)


def test_composite_demand_composition():
    original = WorkspaceSnapshot(demand_scale=1.2)
    snapshot = original.model_dump(mode="json")

    traffic, n1, effective = (
        prepare_traffic_n1_investigation(
            plan=TrafficN1InvestigationPlan(
                traffic_multiplier=1.3
            ),
            workspace=original,
        )
    )

    assert traffic.kind == "traffic_growth"
    assert traffic.traffic_multiplier == 1.3
    assert n1.kind == "n1_resilience"
    assert n1.execution_mode == "fast"
    assert effective.demand_scale == pytest.approx(1.56)
    assert original.model_dump(mode="json") == snapshot


def test_invalid_composite_multiplier():
    with pytest.raises(ValidationError):
        TrafficN1InvestigationPlan(
            traffic_multiplier=-1.0
        )


def test_effective_scale_limit():
    with pytest.raises(ValueError, match="out of range"):
        prepare_traffic_n1_investigation(
            plan=TrafficN1InvestigationPlan(
                traffic_multiplier=2.0
            ),
            workspace=WorkspaceSnapshot(
                demand_scale=600.0
            ),
        )
