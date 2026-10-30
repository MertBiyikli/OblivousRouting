
from __future__ import annotations

import math
from pydantic import BaseModel, ConfigDict, Field, model_validator


class CapacityModification(BaseModel):
    model_config = ConfigDict(extra="forbid")

    edge_id: int = Field(ge=0)
    factor: float = Field(ge=1.0)

    @model_validator(mode="after")
    def validate_factor(self):
        if not math.isfinite(self.factor):
            raise ValueError("Capacity factor must be finite")
        return self


class CommodityModification(BaseModel):
    model_config = ConfigDict(extra="forbid")

    source: int = Field(ge=0)
    target: int = Field(ge=0)
    factor: float = Field(gt=0)

    @model_validator(mode="after")
    def validate_values(self):
        if self.source == self.target:
            raise ValueError("Commodity endpoints must differ")

        if not math.isfinite(self.factor):
            raise ValueError("Commodity factor must be finite")

        return self


class WorkspaceSnapshot(BaseModel):
    model_config = ConfigDict(extra="forbid")

    capacity_upgrades: list[CapacityModification] = Field(
        default_factory=list
    )

    demand_scale: float = Field(default=1.0, gt=0)

    commodity_overrides: list[CommodityModification] = Field(
        default_factory=list
    )

    @model_validator(mode="after")
    def validate_workspace(self):
        if not math.isfinite(self.demand_scale):
            raise ValueError("Demand scale must be finite")

        capacity_ids = [item.edge_id for item in self.capacity_upgrades]

        if len(capacity_ids) != len(set(capacity_ids)):
            raise ValueError("Duplicate capacity upgrade")

        commodity_pairs = [
            (item.source, item.target)
            for item in self.commodity_overrides
        ]

        if len(commodity_pairs) != len(set(commodity_pairs)):
            raise ValueError("Duplicate commodity override")

        return self

    def normalized(self) -> "WorkspaceSnapshot":
        return WorkspaceSnapshot(
            capacity_upgrades=sorted(
                self.capacity_upgrades,
                key=lambda item: item.edge_id,
            ),
            demand_scale=self.demand_scale,
            commodity_overrides=sorted(
                self.commodity_overrides,
                key=lambda item: (item.source, item.target),
            ),
        )
