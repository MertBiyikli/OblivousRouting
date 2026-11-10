
"""Validated experiment specifications for Jeylan."""

from __future__ import annotations

from typing import Literal

from pydantic import BaseModel, ConfigDict, Field


class CapacityExperimentPlan(BaseModel):
    model_config = ConfigDict(extra="forbid")

    edge_id: int = Field(strict=True, ge=0)

    capacity_factor: float = Field(
        gt=1.0,
        le=100.0,
        allow_inf_nan=False,
    )

    demand_model: Literal[
        "uniform",
        "gravity",
        "bimodal",
        "gaussian",
    ] = "uniform"


class CapacityExperimentResult(BaseModel):
    model_config = ConfigDict(extra="forbid")

    network_id: str
    analysis_id: str
    analysis_type: Literal["capacity_upgrade"] = "capacity_upgrade"
    status: Literal["completed"] = "completed"

    edge_id: int
    capacity_factor: float
    demand_model: str

    baseline_max_congestion: float
    fixed_routing_max_congestion: float
    reoptimized_max_congestion: float

    fixed_routing_change: float | None = None
    reoptimized_change: float | None = None

    evidence_source: Literal[
        "persisted_analysis_report"
    ] = "persisted_analysis_report"
