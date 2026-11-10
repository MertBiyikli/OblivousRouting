
"""Validated requests for bounded composite investigations."""

from __future__ import annotations

import math
from typing import Annotated, Literal, Union

from pydantic import (
    BaseModel,
    ConfigDict,
    Field,
    model_validator,
)

from e_routing.api.workspace_models import WorkspaceSnapshot
from e_routing.assistant.analysis_plan import (
    CapacityAnalysisPlan,
    CapacityUpgrade,
    LinkFailureAnalysisPlan,
    N1AnalysisPlan,
    TrafficGrowthAnalysisPlan,
)


class TrafficN1InvestigationPlan(BaseModel):
    model_config = ConfigDict(extra="forbid")

    kind: Literal["traffic_n1"] = "traffic_n1"

    traffic_multiplier: float = Field(
        gt=0.0,
        le=1000.0,
        allow_inf_nan=False,
    )

    demand_model: Literal[
        "uniform", "gravity", "bimodal", "gaussian"
    ] = "uniform"

    execution_mode: Literal[
        "fast", "standard"
    ] = "fast"

    seed: int = 42

    @model_validator(mode="after")
    def validate_multiplier(self):
        if not math.isfinite(self.traffic_multiplier):
            raise ValueError("Non-finite traffic multiplier")
        return self


def prepare_traffic_n1_investigation(
        *,
        plan: TrafficN1InvestigationPlan,
        workspace: WorkspaceSnapshot,
) -> tuple[
    TrafficGrowthAnalysisPlan,
    N1AnalysisPlan,
    WorkspaceSnapshot,
]:
    """Prepare two analyses sharing one effective workspace."""

    original = WorkspaceSnapshot.model_validate(
        workspace.model_dump(mode="json")
    ).normalized()

    effective_scale = (
            original.demand_scale * plan.traffic_multiplier
    )

    if (
            not math.isfinite(effective_scale)
            or effective_scale <= 0
            or effective_scale > 1000
    ):
        raise ValueError(
            "Effective demand multiplier is out of range"
        )

    effective_data = original.model_dump(mode="json")
    effective_data["demand_scale"] = effective_scale

    effective_workspace = WorkspaceSnapshot.model_validate(
        effective_data
    ).normalized()

    traffic = TrafficGrowthAnalysisPlan(
        demand_model=plan.demand_model,
        traffic_multiplier=plan.traffic_multiplier,
    )

    n1 = N1AnalysisPlan(
        demand_model=plan.demand_model,
        execution_mode=plan.execution_mode,
        seed=plan.seed,
    )

    return traffic, n1, effective_workspace



# ------------------------------------------------------------
# Generalized sequential composite analysis
# ------------------------------------------------------------

DemandModel = Literal[
    "uniform", "gravity", "bimodal", "gaussian"
]

ExecutionMode = Literal["fast", "standard"]


# ------------------------------------------------------------
# Generalized composite plans — F3-A
# ------------------------------------------------------------

class CompositeTrafficStep(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=True)

    operation: Literal["traffic_growth"]
    traffic_multiplier: float = Field(
        gt=0,
        le=1000,
        allow_inf_nan=False,
    )


class CompositeCapacityStep(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=True)

    operation: Literal["capacity_upgrade"]
    upgrades: list[CapacityUpgrade] = Field(
        min_length=1,
        max_length=1,
    )


class CompositeN1Step(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=True)

    operation: Literal["n1_resilience"]
    execution_mode: Literal["fast", "standard"] = "fast"


class CompositeLinkFailureStep(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=True)

    operation: Literal["link_failure"]
    failed_edge_id: int = Field(ge=0)


CompositeStep = Annotated[
    CompositeTrafficStep
    | CompositeCapacityStep
    | CompositeN1Step
    | CompositeLinkFailureStep,
    Field(discriminator="operation"),
]


class CompositeAnalysisPlan(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=True)

    kind: Literal["composite"] = "composite"

    demand_model: Literal[
        "uniform", "gravity", "bimodal", "gaussian"
    ] = "uniform"

    # Current Traffic Growth and composite provenance
    # rely on the seeded demand configuration.
    seed: Literal[42] = 42

    steps: list[CompositeStep] = Field(
        min_length=2,
        max_length=6,
    )

    @model_validator(mode="after")
    def validate_sequence(self):
        traffic_count = sum(
            isinstance(step, CompositeTrafficStep)
            for step in self.steps
        )

        capacity_count = sum(
            isinstance(step, CompositeCapacityStep)
            for step in self.steps
        )

        if traffic_count > 1:
            raise ValueError(
                "Only one traffic transition is currently supported"
            )

        if capacity_count > 1:
            raise ValueError(
                "Only one capacity comparison is currently supported"
            )

        return self


def single_analysis_plan(
        step: CompositeStep,
        composite: CompositeAnalysisPlan,
):
    """Translate one composite step to a trusted analysis plan."""

    common = {
        "demand_model": composite.demand_model,
    }

    if isinstance(step, CompositeTrafficStep):
        return TrafficGrowthAnalysisPlan(
            **common,
            traffic_multiplier=step.traffic_multiplier,
        )

    if isinstance(step, CompositeCapacityStep):
        return CapacityAnalysisPlan(
            **common,
            upgrades=step.upgrades,
        )

    if isinstance(step, CompositeN1Step):
        return N1AnalysisPlan(
            **common,
            execution_mode=step.execution_mode,
            seed=composite.seed,
        )

    if isinstance(step, CompositeLinkFailureStep):
        return LinkFailureAnalysisPlan(
            **common,
            failed_edge_id=step.failed_edge_id,
        )

    raise ValueError("Unsupported composite operation")
