
"""Network-scoped fixed-routing Traffic Growth API."""

from __future__ import annotations

import json
import math
import os
import re
import subprocess
from pathlib import Path
from typing import Literal

from fastapi import APIRouter, HTTPException
from fastapi.concurrency import run_in_threadpool
from pydantic import BaseModel, ConfigDict, Field, model_validator
import logging

logger = logging.getLogger(__name__)
from e_routing.api.config import (
    load_network_registry,
    trusted_run_for_network,
)

_ID = re.compile(r"^[A-Za-z0-9_-]{1,64}$")


class CommodityOverride(BaseModel):
    model_config = ConfigDict(extra="forbid")

    source: int = Field(ge=0)
    target: int = Field(ge=0)
    factor: float = Field(gt=0, le=1000, allow_inf_nan=False)

    @model_validator(mode="after")
    def validate_pair(self):
        if self.source == self.target:
            raise ValueError("Source and target must differ")
        return self


class TrafficRequest(BaseModel):
    model_config = ConfigDict(extra="forbid")

    analysis_id: str = Field(pattern=r"^[A-Za-z0-9_-]{1,64}$")
    demand_model: Literal["uniform"] = "uniform"
    seed: int = Field(default=42, ge=0, le=2147483647)
    global_scale: float = Field(
        default=1.0, gt=0, le=1000, allow_inf_nan=False
    )
    commodity_overrides: list[CommodityOverride] = Field(
        default_factory=list, max_length=1000
    )

    @model_validator(mode="after")
    def unique_commodities(self):
        pairs = [(x.source, x.target) for x in self.commodity_overrides]
        if len(pairs) != len(set(pairs)):
            raise ValueError("Duplicate commodity overrides")
        return self


def validate_traffic_report(report: object) -> dict:
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
        raise ValueError("Baseline and modified links differ")

    return report


def create_traffic_router(
        reports_root: Path,
        timeout_seconds: int = 120,
) -> APIRouter:
    router = APIRouter()
    root = Path(reports_root).resolve()


    def write_status(
            folder: Path,
            analysis_id: str,
            state: str,
            message: str | None = None,
    ) -> None:
        payload = {
            "analysis_id": analysis_id,
            "status": state,
        }
        if message is not None:
            payload["message"] = message

        temporary = folder / "status.pending.json"
        destination = folder / "status.json"

        temporary.write_text(
            json.dumps(payload),
            encoding="utf-8",
        )
        os.replace(temporary, destination)




    def record_failure(
            folder: Path,
            analysis_id: str,
            message: str,
    ) -> None:
        try:
            write_status(
                folder,
                analysis_id,
                "failed",
                message,
            )
        except OSError:
            logger.exception(
                "Could not persist failed status for %s",
                analysis_id,
            )



    def resolve_network(network_id: str):
        if not _ID.fullmatch(network_id):
            raise HTTPException(422, "Invalid network ID")

        try:
            registry = load_network_registry()
        except (OSError, ValueError, RuntimeError):
            raise HTTPException(503, "Network registry unavailable")

        network = registry.get(network_id)
        if network is None:
            raise HTTPException(404, "Network not found")

        try:
            trusted = trusted_run_for_network(network)
        except (OSError, ValueError, RuntimeError):
            raise HTTPException(503, "Network artifacts unavailable")

        return trusted

    def report_paths(network_id: str, analysis_id: str):
        if not _ID.fullmatch(analysis_id):
            raise HTTPException(422, "Invalid analysis ID")

        folder = root / network_id / analysis_id
        return (
            folder,
            folder / "traffic.json",
            folder / "traffic.pending.json",
        )

    @router.post(
        "/api/v1/networks/{network_id}/analyses/traffic",
        status_code=201,
    )
    async def create_traffic(
            network_id: str,
            request: TrafficRequest,
    ):
        trusted = await run_in_threadpool(
            resolve_network, network_id
        )

        folder, report_path, pending_path = report_paths(
            network_id, request.analysis_id
        )
        folder.mkdir(parents=True, exist_ok=True)

        reservation = folder / "traffic_request.json"

        try:
            fd = os.open(
                reservation,
                os.O_WRONLY | os.O_CREAT | os.O_EXCL,
                0o600,
                )
        except FileExistsError:
            raise HTTPException(409, "Analysis ID already exists")
        except OSError:
            raise HTTPException(503, "Analysis storage unavailable")

        with os.fdopen(fd, "w", encoding="utf-8") as output:
            json.dump(
                {
                    "network_id": network_id,
                    **request.model_dump(),
                },
                output,
            )

        try:
            write_status(
                    folder,
                    request.analysis_id,
                    "running",
            )
        except OSError:
            logger.exception("Unable to persist traffic analysis state")
            raise HTTPException(503, "Analysis storage unavailable")

        command = [
            str(trusted.executable),
            "analyze-traffic",
            "--solver", trusted.solver_label,
            "--graph", str(trusted.graph),
            "--load-routing-table", str(trusted.routing_table),
            "--demand", request.demand_model,
            "--seed", str(request.seed),
            "--demand-scale", str(request.global_scale),
        ]

        if request.commodity_overrides:
            command.append("--commodity-scale")
            command.extend(
                f"{x.source}:{x.target}:{x.factor}"
                for x in request.commodity_overrides
            )

        command.extend(["--output", str(pending_path)])

        def execute():
            return subprocess.run(
                command,
                capture_output=True,
                text=True,
                timeout=timeout_seconds,
                check=False,
            )


        try:
            completed = await run_in_threadpool(execute)

        except subprocess.TimeoutExpired:
            logger.exception("Traffic analysis timed out")
            record_failure(
                folder,
                request.analysis_id,
                "Traffic engine returned an error",
            )
            raise HTTPException(504, "Traffic analysis timed out")

        except OSError:
            logger.exception("Unable to execute traffic engine")
            record_failure(
                folder,
                request.analysis_id,
                "Traffic engine returned an error",
            )
            raise HTTPException(503, "Traffic engine unavailable")


        if completed.returncode != 0:
            logger.error(
                "Traffic analysis %s failed (exit=%s)\n"
                "Command: %s\n"
                "stdout: %s\n"
                "stderr: %s",
                request.analysis_id,
                completed.returncode,
                command,
                completed.stdout[-4000:],
                completed.stderr[-4000:],
            )
            record_failure(
                folder,
                request.analysis_id,
                "Traffic engine returned an error",
            )
            raise HTTPException(502, "Traffic analysis failed")

        try:
            report = validate_traffic_report(
                json.loads(pending_path.read_text(encoding="utf-8"))
            )

            if not math.isclose(
                    report["requested_global_scale"],
                    request.global_scale,
                    rel_tol=1e-10,
            ):
                raise ValueError("Global scale mismatch")

            expected = sorted(
                (
                    x.source,
                    x.target,
                    x.factor,
                )
                for x in request.commodity_overrides
            )

            actual = sorted(
                (
                    x["source"],
                    x["target"],
                    x["factor"],
                )
                for x in report["commodity_overrides"]
            )

            if actual != expected:
                raise ValueError("Commodity overrides mismatch")

            pending_path.replace(report_path)

        except (OSError, ValueError) as exc:
            logger.exception("Traffic analysis report validation failed")
            record_failure(
                folder,
                request.analysis_id,
                "Traffic engine returned an error",
            )
            raise HTTPException(
            502, "Invalid traffic analysis report"
            ) from exc

        try:
            write_status(
                folder,
                request.analysis_id,
                "completed",
            )
        except OSError:
            logger.exception("Unable to persist completed state")
            raise HTTPException(503, "Analysis storage unavailable")

        return {
            "network_id": network_id,
            "analysis_id": request.analysis_id,
            "status": "completed",
            "execution_mode": "fixed_routing",
            "summary": {
                "baseline_max_congestion":
                    report["baseline"]["max_congestion"],
                "modified_max_congestion":
                    report["modified"]["max_congestion"],
                "baseline_demand_volume":
                    report["baseline"]["total_demand_volume"],
                "modified_demand_volume":
                    report["modified"]["total_demand_volume"],
            },
        }

    @router.get(
        "/api/v1/networks/{network_id}/analyses/traffic/{analysis_id}"
    )
    async def get_traffic(network_id: str, analysis_id: str):
        await run_in_threadpool(resolve_network, network_id)



        folder, report_path, _ = report_paths(
            network_id, analysis_id
        )

        reservation = folder / "traffic_request.json"
        status_path = folder / "status.json"

        if not reservation.is_file():
            raise HTTPException(404, "Traffic analysis not found")

        if not report_path.is_file():
            if status_path.is_file():
                try:
                    state = json.loads(
                        status_path.read_text(encoding="utf-8")
                    )
                except (OSError, ValueError):
                    raise HTTPException(
                        500, "Analysis status invalid"
                    )

                if state.get("analysis_id") != analysis_id:
                    raise HTTPException(
                        500, "Analysis ID mismatch"
                    )

                if state.get("status") in {"running", "failed"}:
                    return state

            return {
                "analysis_id": analysis_id,
                "status": "running",
            }


        try:
            return validate_traffic_report(
                json.loads(report_path.read_text(encoding="utf-8"))
            )
        except (OSError, ValueError, json.JSONDecodeError):
            raise HTTPException(500, "Stored traffic report invalid")

    return router
