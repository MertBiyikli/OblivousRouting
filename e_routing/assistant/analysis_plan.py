
"""Validated analysis plans shared by Jeylan execution workflows."""

from __future__ import annotations

from typing import Annotated, Literal

from pydantic import (
    BaseModel,
    ConfigDict,
    Field,
    model_validator,
)


DemandModel = Literal[
    "uniform",
    "gravity",
    "bimodal",
    "gaussian",
]


class AnalysisPlan(BaseModel):
    model_config = ConfigDict(
        extra="forbid",
        frozen=True,
    )

    demand_model: DemandModel = "uniform"


class CapacityUpgrade(BaseModel):
    model_config = ConfigDict(
        extra="forbid",
        frozen=True,
    )

    edge_id: int = Field(ge=0)
    capacity_factor: float = Field(
        gt=1.0,
        le=100.0,
        allow_inf_nan=False,
    )


class CapacityAnalysisPlan(AnalysisPlan):
    kind: Literal["capacity_upgrade"] = "capacity_upgrade"

    # The current detailed-analysis workflow supports one
    # event upgrade for a three-way capacity comparison.
    upgrades: list[CapacityUpgrade] = Field(
        min_length=1,
        max_length=1,
    )


class LinkFailureAnalysisPlan(AnalysisPlan):
    kind: Literal["link_failure"] = "link_failure"

    failed_edge_id: int = Field(ge=0)


class N1AnalysisPlan(AnalysisPlan):
    kind: Literal["n1_resilience"] = "n1_resilience"

    execution_mode: Literal[
        "fast",
        "standard",
    ] = "fast"

    seed: int = Field(
        default=42,
        ge=0,
        le=2147483647,
    )


class TrafficGrowthAnalysisPlan(AnalysisPlan):
    kind: Literal["traffic_growth"] = "traffic_growth"

    # Relative multiplier applied to the existing workspace
    # demand configuration.
    traffic_multiplier: float = Field(
        gt=0,
        le=1000,
        allow_inf_nan=False,
    )


AnalysisOperation = Annotated[
    CapacityAnalysisPlan
    | LinkFailureAnalysisPlan
    | N1AnalysisPlan
    | TrafficGrowthAnalysisPlan,
    Field(discriminator="kind"),
]
