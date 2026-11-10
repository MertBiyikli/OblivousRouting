"""Opt-in C++ scenario recomputation; mount alongside existing fast endpoints.
Requires the C++ `ER_SCENARIO_REPORT_PATH` hook described in INTEGRATION.md.
"""
from __future__ import annotations

import json
import logging
import os
import re
import subprocess
from pathlib import Path
from typing import Literal

from fastapi import APIRouter, HTTPException
from pydantic import BaseModel, ConfigDict, Field, model_validator

from e_routing.api.config import load_network_registry
from e_routing.api.scenario_report_api import build_analysis_payload
from e_routing.api.workspace_cli import workspace_cli_args
from e_routing.api.workspace_models import WorkspaceSnapshot

_ID = re.compile(r"^[A-Za-z0-9_-]{1,80}$")

logger = logging.getLogger(__name__)


class Upgrade(BaseModel):
    edge_id: int = Field(ge=0)
    capacity_factor: float = Field(
        ge=1.0, le=100.0, allow_inf_nan=False
    )


from typing import Literal
from pydantic import BaseModel, Field, ConfigDict, model_validator
from e_routing.api.workspace_models import WorkspaceSnapshot


class DetailedRequest(BaseModel):
    model_config = ConfigDict(extra="forbid")

    scenario_id: str = Field(
        pattern=r"^[A-Za-z0-9_-]{1,64}$"
    )

    scenario_kind: Literal[
        "capacity_upgrade", "link_failure"
    ] = "capacity_upgrade"

    demand_model: Literal[
        "uniform", "gravity", "bimodal", "gaussian"
    ] = "uniform"

    workspace: WorkspaceSnapshot = Field(
        default_factory=WorkspaceSnapshot
    )

    # Event upgrades, not pre-existing workspace upgrades.
    capacity_upgrades: list[Upgrade] = Field(
        default_factory=list,
        max_length=100
    )

    link_failures: list[int] = Field(
        default_factory=list,
        max_length=1
    )

    @model_validator(mode="after")
    def validate_modification(self):
        if self.scenario_kind == "link_failure":
            if len(self.link_failures) != 1:
                raise ValueError(
                    "Select exactly one failed physical link"
                )

            if self.link_failures[0] < 0:
                raise ValueError(
                    "Failed edge ID must be nonnegative"
                )

            if self.capacity_upgrades:
                raise ValueError(
                    "Link-failure events cannot include capacity upgrades"
                )

        elif self.scenario_kind == "capacity_upgrade":
            if self.link_failures:
                raise ValueError(
                    "Capacity analysis cannot fail links"
                )

            if not self.capacity_upgrades:
                raise ValueError(
                    "Select at least one capacity upgrade"
                )

            # Current capacity CLI supports only one event upgrade.
            # It cannot also apply pre-existing workspace upgrades.
            if self.workspace.capacity_upgrades:
                raise ValueError(
                    "Capacity comparison with existing workspace "
                    "capacity upgrades is not yet supported"
                )

        return self


def detailed_report_path(root: Path, network_id: str, scenario_id: str) -> Path:
    if not _ID.fullmatch(network_id) or not _ID.fullmatch(scenario_id):
        raise HTTPException(422, "Invalid network/scenario ID")
    return root / "detailed" / network_id / scenario_id / "report.json"



def execute_detailed_analysis(
        network_id: str,
        request: DetailedRequest,
        reports_root: Path | str,
) -> dict:
    root = Path(reports_root).resolve()
    path = detailed_report_path(root, network_id, request.scenario_id)

    registry = load_network_registry()
    network = registry.get(network_id)
    if network is None:
        raise HTTPException(404, "Unknown registered network")

    network.verify()

    if request.demand_model not in {
        "uniform", "gravity", "bimodal", "gaussian"
    }:
        raise HTTPException(422, "Unsupported demand model")

    ids = [u.edge_id for u in request.capacity_upgrades]
    if len(ids) != len(set(ids)):
        raise HTTPException(422, "Duplicate upgrade edge IDs")

    # Server-controlled output paths.
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() or (path.parent / "request.json").exists():
        raise HTTPException(409, "Scenario ID already exists")

    executable = os.environ.get("ER_ENGINE")
    if not executable:
        raise HTTPException(503, "ER_ENGINE is not configured")

    engine = Path(executable).resolve()
    if not engine.is_file():
        raise HTTPException(503, "Configured C++ engine not found")

    metadata = {
        "network_id": network_id,
        "scenario_id": request.scenario_id,
        "scenario_kind": request.scenario_kind,
        "demand_model": request.demand_model,
        "workspace": request.workspace.model_dump(mode="json"),
        "capacity_upgrades": [
            u.model_dump() for u in request.capacity_upgrades
        ],
        "link_failures": request.link_failures,
    }

    # Exclusive reservation prevents scenario ID reuse.
    try:
        with (path.parent / "request.json").open(
                "x", encoding="utf-8"
        ) as handle:
            json.dump(metadata, handle)
    except FileExistsError:
        raise HTTPException(409, "Scenario ID already exists")

    analysis_kind = (
        "capacity"
        if request.scenario_kind == "capacity_upgrade"
        else "link-failure"
    )

    command = [
        str(engine),
        "analyze-scenario",
        "--analysis", analysis_kind,
        "--solver", network.solver_label,
        "--graph", str(network.graph),
        "--demand", request.demand_model,
        "--seed", "42",
        "--execution-mode", "fast",
    ]

    command.extend(
        workspace_cli_args(
            request.workspace,
            include_capacity=(
                    request.scenario_kind == "link_failure"
            ),
        )
    )

    if request.scenario_kind == "capacity_upgrade":
        for upgrade in request.capacity_upgrades:
            command.extend([
                "--capacity-upgrade",
                f"{upgrade.edge_id}:{upgrade.capacity_factor}",
            ])
    else:
        command.extend([
            "--link-failure",
            str(request.link_failures[0]),
        ])

    command.extend(["--output", str(path)])

    try:
        completed = subprocess.run(
            command,
            env=os.environ.copy(),
            cwd=str(root.parent.parent),
            capture_output=True,
            text=True,
            timeout=180,
        )
    except subprocess.TimeoutExpired:
        raise HTTPException(
            504,
            "Detailed C++ scenario exceeded timeout",
        )
    except OSError as exc:
        raise HTTPException(
            502,
            f"Could not launch C++ engine: {exc}",
        ) from exc

    if completed.returncode != 0:
        stderr = completed.stderr or ""

        logger.error(
            "Detailed scenario %s failed (exit=%s)\n"
            "stdout: %s\n"
            "stderr: %s",
            request.scenario_id,
            completed.returncode,
            completed.stdout[-4000:],
            stderr[-4000:],
        )

        raise HTTPException(
            502,
            "C++ detailed scenario execution failed",
        )

    if not path.is_file():
        raise HTTPException(
            502,
            "C++ did not write the requested detailed report",
        )

    try:
        report = json.loads(path.read_text(encoding="utf-8"))

        if not isinstance(report, dict):
            raise ValueError("Expected a JSON object")

        solvers = report.get("solver_results")

        is_disconnected = (
                report.get("schema_version") == "1.3"
                and report.get("scenario_kind") == "link_failure"
                and isinstance(solvers, list)
                and len(solvers) > 0
                and all(
            isinstance(solver, dict)
            and isinstance(solver.get("baseline"), dict)
            and solver.get("reoptimized") is None
            and solver.get("comparison") is None
            for solver in solvers
        )
        )

        if is_disconnected:
            if request.scenario_kind != "link_failure":
                raise ValueError(
                    "Unexpected disconnected analysis outcome"
                )

            return {
                "scenario_id": request.scenario_id,
                "network_id": network_id,
                "status": "disconnected",
                "report_available": True,
                "failed_edge_id": request.link_failures[0],
                "reason": "Link failure disconnects the network",
            }

        if report.get("status", "completed") != "completed":
            raise ValueError(
                "Unknown scenario execution status"
            )

        payload = build_analysis_payload(report)

    except (ValueError, KeyError, TypeError) as exc:
        raise HTTPException(
            502,
            f"Invalid C++ detailed report: {exc}",
        ) from exc

    return {
        "scenario_id": request.scenario_id,
        "network_id": network_id,
        "status": "completed",
        "report_available": True,
        "analysis": payload,
    }

def create_detailed_analysis_router(
        reports_root: Path | str = "results/scenarios",
) -> APIRouter:
    router = APIRouter()
    root = Path(reports_root).resolve()

    @router.post(
        "/api/v1/networks/{network_id}/scenarios/detailed-analysis"
    )
    def run_detailed_analysis(
            network_id: str,
            request: DetailedRequest,
    ):
        return execute_detailed_analysis(
            network_id=network_id,
            request=request,
            reports_root=root,
        )

    return router
