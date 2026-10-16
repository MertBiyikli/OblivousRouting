
"""HTTP interface for deterministic E-Routing scenario execution."""

from __future__ import annotations

import logging
import os
from dataclasses import replace
from pathlib import Path
from typing import Any
from uuid import uuid4

from fastapi import FastAPI, HTTPException
from fastapi.concurrency import run_in_threadpool
from fastapi.responses import FileResponse, RedirectResponse

from e_routing.api.config import load_trusted_config
from e_routing.scenarios.repository import ScenarioRepository
from e_routing.scenarios.runtime import ScenarioRequest, execute

logger = logging.getLogger(__name__)

app = FastAPI(title="E-Routing API", version="0.1.0")


@app.get("/", include_in_schema=False)
def home() -> RedirectResponse:
    return RedirectResponse(url="/ui")

@app.get("/ui", include_in_schema=False)
def web_ui() -> FileResponse:
    return FileResponse(Path(__file__).parent / "static" / "index.html")

@app.get("/health")
def health() -> dict[str, str]:
    return {"status": "ok", "service": "e-routing"}


@app.post("/api/v1/scenarios")

@app.post("/api/v1/scenarios")
async def create_scenario(payload: dict[str, Any]) -> dict[str, Any]:
    try:
        request = ScenarioRequest.from_dict(payload)
    except ValueError as exc:
        raise HTTPException(status_code=422, detail=str(exc)) from exc

    try:
        trusted = load_trusted_config()
        repository = ScenarioRepository(Path(os.environ.get("ER_RESPONSE_DIR", "results/scenario_responses")))
    except (RuntimeError, OSError) as exc:
        logger.exception("Scenario API is not configured")
        raise HTTPException(status_code=503, detail="Scenario service unavailable") from exc

    try:
        await run_in_threadpool(repository.reserve, request.scenario_id)
    except FileExistsError as exc:
        raise HTTPException(status_code=409, detail="Scenario ID already exists") from exc
    except (OSError, ValueError) as exc:
        logger.exception("Scenario reservation failed")
        raise HTTPException(status_code=503, detail="Scenario storage unavailable") from exc

    execution_id = uuid4().hex
    isolated_trusted = replace(trusted, workspace=trusted.workspace / execution_id)

    try:
        result = await run_in_threadpool(execute, request, isolated_trusted)
        await run_in_threadpool(repository.complete, result)
        return result
    except (OSError, RuntimeError, ValueError) as exc:
        logger.exception("Scenario execution failed")
        try:
            await run_in_threadpool(repository.release, request.scenario_id)
        except Exception:
            logger.exception("Could not release failed scenario reservation")
        raise HTTPException(status_code=502, detail="Scenario execution failed") from exc




@app.get("/api/v1/scenarios/{scenario_id}")
def get_scenario(scenario_id: str) -> dict[str, Any]:
    repository = ScenarioRepository(
        Path(os.environ.get("ER_RESPONSE_DIR", "results/scenario_responses"))
    )

    try:
        result = repository.get(scenario_id)
    except ValueError as exc:
        raise HTTPException(status_code=422, detail=str(exc)) from exc

    if result is None:
        raise HTTPException(status_code=404, detail="Scenario not found")

    if result.get("status") == "running":
        raise HTTPException(status_code=409, detail="Scenario is still running")

    return result
