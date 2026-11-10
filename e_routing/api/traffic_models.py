
"""Traffic Growth request models and report validation."""

from __future__ import annotations

import math
from typing import Literal

from pydantic import (
    BaseModel,
    ConfigDict,
    Field,
    model_validator,
)

from e_routing.api.detailed_analysis_api import Upgrade


class CommodityOverride(BaseModel):
    model_config = ConfigDict(extra="forbid")

    source: int = Field(ge=0)
    target: int = Field(ge=0)
    factor: float = Field(
        gt=0,
        le=1000,
        allow_inf_nan=False,
    )

    @model_validator(mode="after")
    def validate_pair(self):
        if self.source == self.target:
            raise ValueError("Source and target must differ")
        return self


class TrafficRequest(BaseModel):
    model_config = ConfigDict(extra="forbid")

    analysis_id: str = Field(
        pattern=r"^[A-Za-z0-9_-]{1,64}$"
    )
    demand_model: Literal[
        "uniform",
        "gravity",
        "bimodal",
        "gaussian",
    ] = "uniform"
    seed: int = Field(
        default=42,
        ge=0,
        le=2147483647,
    )
    global_scale: float = Field(
        default=1.0,
        gt=0,
        le=1000,
        allow_inf_nan=False,
    )
    commodity_overrides: list[CommodityOverride] = Field(
        default_factory=list,
        max_length=1000,
    )
    capacity_upgrades: list[Upgrade] = Field(
        default_factory=list,
        max_length=100,
    )

    @model_validator(mode="after")
    def unique_commodities(self):
        pairs = [
            (item.source, item.target)
            for item in self.commodity_overrides
        ]
        ids = [
            upgrade.edge_id
            for upgrade in self.capacity_upgrades
        ]

        if len(ids) != len(set(ids)):
            raise ValueError(
                "Duplicate capacity upgrade edge IDs"
            )

        if len(pairs) != len(set(pairs)):
            raise ValueError("Duplicate commodity overrides")

        return self


def validate_traffic_report(report: object) -> dict:
    """Validate the persisted fixed-routing Traffic Growth report."""

    if not isinstance(report, dict):
        raise ValueError("Report must be a JSON object")

    if report.get("scenario_kind") != "traffic_growth":
        raise ValueError("Incorrect scenario kind")

    if report.get("execution_mode") != "fixed_routing":
        raise ValueError("Incorrect execution mode")

    def number(value):
        return (
                type(value) in (int, float)
                and math.isfinite(value)
        )

    for field in (
            "demand_scale",
            "requested_global_scale",
            "absolute_demand_change",
            "absolute_congestion_change",
    ):
        if not number(report.get(field)):
            raise ValueError(f"Invalid metric: {field}")

    for field in (
            "relative_demand_change",
            "relative_congestion_change",
    ):
        value = report.get(field)
        if value is not None and not number(value):
            raise ValueError(f"Invalid metric: {field}")

    overrides = report.get("commodity_overrides")
    if not isinstance(overrides, list):
        raise ValueError("Missing commodity overrides")

    for item in overrides:
        if not isinstance(item, dict):
            raise ValueError("Invalid commodity override")

        if (
                type(item.get("source")) is not int
                or type(item.get("target")) is not int
                or not number(item.get("factor"))
                or item["factor"] <= 0
        ):
            raise ValueError("Invalid commodity override values")

    link_sets = []

    for name in ("baseline", "modified"):
        network = report.get(name)

        if not isinstance(network, dict):
            raise ValueError(f"Missing {name} network")

        for metric in (
                "total_demand_volume",
                "max_congestion",
                "average_congestion",
        ):
            if not number(network.get(metric)):
                raise ValueError(f"Invalid {name}.{metric}")

        links = network.get("links")
        if not isinstance(links, list) or not links:
            raise ValueError(f"Missing {name} links")

        ids = set()

        for link in links:
            if not isinstance(link, dict):
                raise ValueError("Invalid link")

            edge_id = link.get("edge_id")

            if type(edge_id) is not int or edge_id in ids:
                raise ValueError("Invalid or duplicate edge ID")

            ids.add(edge_id)

            for metric in ("capacity", "load", "utilization"):
                value = link.get(metric)
                if not number(value) or value < 0:
                    raise ValueError(f"Invalid link {metric}")

            if link["capacity"] <= 0:
                raise ValueError("Nonpositive link capacity")

        link_sets.append(ids)

    if link_sets[0] != link_sets[1]:
        raise ValueError(
            "Baseline and modified links differ"
        )

    return report
