
"""Network-scoped N-1 resilience API."""

from __future__ import annotations

import asyncio
import json
import logging
import os
import re
from pathlib import Path
from typing import Literal

from fastapi import APIRouter, HTTPException
from fastapi.concurrency import run_in_threadpool
from pydantic import BaseModel, ConfigDict, Field
from e_routing.api.workspace_models import WorkspaceSnapshot
from e_routing.api.config import (
    load_network_registry,
    trusted_run_for_network,
)
from e_routing.api.n1_report_api import validate_n1_report

logger = logging.getLogger(__name__)

_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_-]{0,63}$")

# Explicit CLI-to-report solver identities.
_SOLVERS = {
    "electrical": "electrical_sketching",
}


class N1Request(BaseModel):
    model_config = ConfigDict(extra="forbid")

    analysis_id: str = Field(min_length=1, max_length=64)
    demand_model: Literal["uniform"] = "uniform"
    seed: int = Field(default=42, ge=0, le=2147483647)
    execution_mode: Literal[
        "fast", "standard", "research"
    ] = "fast"
    workspace: WorkspaceSnapshot = Field(
        default_factory=WorkspaceSnapshot
    )


def create_n1_router(
        reports_root: Path,
        timeout_seconds: float = 600.0,
) -> APIRouter:

    router = APIRouter()
    root = Path(reports_root).resolve()

    def validate_id(value: str) -> str:
        if not _ID.fullmatch(value):
            raise HTTPException(422, "Invalid analysis or network ID")
        return value

    def resolve_network(network_id: str):
        validate_id(network_id)

        try:
            registry = load_network_registry()
        except (OSError, RuntimeError, ValueError) as exc:
            logger.exception("N-1 registry unavailable")
            raise HTTPException(
                503, "Network registry unavailable"
            ) from exc

        network = registry.get(network_id)

        if network is None:
            raise HTTPException(404, "Network not found")

        if network.solver_label not in _SOLVERS:
            raise HTTPException(
                422, "N-1 solver is not supported"
            )

        try:
            trusted = trusted_run_for_network(network)
            executable = trusted.executable.resolve(strict=True)
        except (OSError, RuntimeError, ValueError) as exc:
            logger.exception("N-1 network verification failed")
            raise HTTPException(
                503, "Network artifacts unavailable"
            ) from exc

        return network, trusted, executable

    def analysis_path(network_id: str, analysis_id: str):
        validate_id(network_id)
        validate_id(analysis_id)
        return root / network_id / analysis_id

    @router.post(
        "/api/v1/networks/{network_id}/analyses/n1",
        status_code=201,
    )
    async def create_n1(network_id: str, request: N1Request):
        if request.execution_mode == "research":
            raise HTTPException(
                status_code=422,
                detail=(
                    "N-1 RESEARCH mode is not available yet. "
                    "Multi-solver analysis must be implemented first."
                ),
            )
        validate_id(request.analysis_id)

        network, trusted, executable = await run_in_threadpool(
            resolve_network, network_id
        )

        normalized_workspace = request.workspace.normalized()



        directory = analysis_path(network_id, request.analysis_id)
        report_path = directory / "n1.json"
        pending_path = directory / "n1.pending.json"
        reservation = directory / "request.json"

        try:
            directory.mkdir(parents=True, exist_ok=True)

            # Atomic reservation across concurrent requests.
            fd = os.open(
                reservation,
                os.O_CREAT | os.O_EXCL | os.O_WRONLY,
                0o600,
                )
            with os.fdopen(fd, "w", encoding="utf-8") as stream:
                json.dump({
                    "network_id": network_id,
                    "analysis_id": request.analysis_id,
                    "status": "running",
                    "request": {
                        **request.model_dump(exclude={"workspace"}),
                        "workspace": normalized_workspace.model_dump(),
                    },
                }, stream)

        except FileExistsError as exc:
            raise HTTPException(
                409, "Analysis ID already exists"
            ) from exc
        except OSError as exc:
            logger.exception("N-1 reservation failed")
            raise HTTPException(
                503, "Analysis storage unavailable"
            ) from exc


        command = [
            str(executable),
            "analyze-n1",
            "--solver", network.solver_label,
            "--graph", str(trusted.graph),
            "--demand", request.demand_model,
            "--seed", str(request.seed),
            "--execution-mode", request.execution_mode,
            "--demand-scale", str(normalized_workspace.demand_scale),
        ]

        for upgrade in normalized_workspace.capacity_upgrades:
            command.extend([
                "--capacity-upgrade",
                f"{upgrade.edge_id}:{upgrade.factor}",
            ])

        if normalized_workspace.commodity_overrides:
            command.append("--commodity-scale")
            command.extend(
                f"{item.source}:{item.target}:{item.factor}"
                for item in normalized_workspace.commodity_overrides
            )

        command.extend([
            "--output", str(pending_path),
        ])


        try:
            process = await asyncio.create_subprocess_exec(
                *command,
                stdout=asyncio.subprocess.PIPE,
                stderr=asyncio.subprocess.PIPE,
            )

            try:
                stdout, stderr = await asyncio.wait_for(
                    process.communicate(),
                    timeout=timeout_seconds,
                )
            except asyncio.TimeoutError:
                process.kill()
                await process.communicate()
                raise HTTPException(
                    504, "N-1 analysis timed out"
                )

            if process.returncode != 0:
                logger.error(
                    "N-1 execution failed: %s",
                    stderr.decode(errors="replace")[-2000:],
                )
                raise HTTPException(
                    502, "N-1 engine execution failed"
                )

            report = validate_n1_report(
                json.loads(pending_path.read_text("utf-8"))
            )

            if report["solver"] != _SOLVERS[network.solver_label]:
                raise HTTPException(
                    502, "N-1 solver identity mismatch"
                )

            # Publish only a validated, complete report.
            os.replace(pending_path, report_path)

            return {
                "network_id": network_id,
                "analysis_id": request.analysis_id,
                "status": "completed",
                "execution_mode": request.execution_mode,
                "summary": report["summary"],
            }

        except HTTPException:
            raise
        except OSError as exc:
            logger.exception("N-1 process or storage failed")
            raise HTTPException(
                503, "N-1 execution unavailable"
            ) from exc
        except (ValueError, KeyError, TypeError) as exc:
            logger.exception("Invalid N-1 report")
            raise HTTPException(
                502, "N-1 report validation failed"
            ) from exc

    @router.get(
        "/api/v1/networks/{network_id}/analyses/n1/{analysis_id}"
    )
    async def get_n1(network_id: str, analysis_id: str):
        await run_in_threadpool(resolve_network, network_id)

        directory = analysis_path(network_id, analysis_id)
        report_path = directory / "n1.json"
        reservation = directory / "request.json"

        if not report_path.is_file():
            if reservation.is_file():
                raise HTTPException(
                    409, "N-1 analysis is not completed"
                )
            raise HTTPException(
                404, "N-1 analysis not found"
            )

        try:
            return await run_in_threadpool(
                lambda: validate_n1_report(
                    json.loads(report_path.read_text("utf-8"))
                )
            )
        except (OSError, ValueError, KeyError, TypeError) as exc:
            logger.exception("Stored N-1 report invalid")
            raise HTTPException(
                500, "Stored N-1 report is invalid"
            ) from exc

    return router
