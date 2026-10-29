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
from pydantic import BaseModel, Field

from e_routing.api.config import load_trusted_config, load_network_registry, trusted_run_for_network
from e_routing.assistant.scenario_chat import answer_scenario_question
from e_routing.scenarios.repository import ScenarioRepository
from e_routing.scenarios.runtime import ScenarioRequest, execute
from e_routing.scenarios.runtime import ScenarioRequest, execute, execute_comparison
from e_routing.api.scenario_report_api import (
    create_scenario_analysis_router,
)

REPORTS_ROOT = Path("results/scenarios").resolve()
logger = logging.getLogger(__name__)

app = FastAPI(title="E-Routing API", version="0.1.0")


from e_routing.api.detailed_analysis_api import (
    create_detailed_analysis_router,
)

app.include_router(
    create_detailed_analysis_router("results/scenarios")
)


import json
from pathlib import Path

REPORTS_ROOT = Path("results/scenarios").resolve()

def lookup_scenario_report(network_id: str, scenario_id: str):
    from e_routing.api.detailed_analysis_api import detailed_report_path

    report = detailed_report_path(
        REPORTS_ROOT, network_id, scenario_id
    )

    request_file = report.parent / "request.json"

    if not report.is_file() or not request_file.is_file():
        return None

    try:
        metadata = json.loads(
            request_file.read_text(encoding="utf-8")
        )
    except (ValueError, OSError):
        return None

    if (
            metadata.get("network_id") != network_id
            or metadata.get("scenario_id") != scenario_id
    ):
        return None

    return report



app.include_router(
    create_scenario_analysis_router(
        reports_root=REPORTS_ROOT,
        lookup_report=lookup_scenario_report,
    )
)


class ScenarioQuestion(BaseModel):
    question: str = Field(min_length=1, max_length=2000)


@app.get("/", include_in_schema=False)
def home() -> RedirectResponse:
    return RedirectResponse(url="/ui")


from fastapi.responses import FileResponse

@app.get("/scenario_studio_analysis.js", include_in_schema=False)
def scenario_studio_analysis_js():
    return FileResponse(
        "e_routing/api/static/scenario_studio_analysis.js",
        media_type="application/javascript"
    )

@app.get("/ui", include_in_schema=False)
def web_ui() -> FileResponse:
    return FileResponse(Path(__file__).parent / "static" / "index.html")


@app.get("/health")
def health() -> dict[str, str]:
    return {"status": "ok", "service": "e-routing"}


@app.post("/api/v1/scenarios")
async def create_scenario(payload: dict[str, Any]) -> dict[str, Any]:
    try:
        ScenarioRequest.from_dict(payload)
    except ValueError as exc:
        raise HTTPException(status_code=422, detail=str(exc)) from exc

    try:
        trusted = load_trusted_config()
    except (RuntimeError, OSError) as exc:
        logger.exception("Scenario API is not configured")
        raise HTTPException(
            status_code=503,
            detail="Scenario service unavailable"
        ) from exc

    return await _run_scenario(payload, trusted)


async def _run_scenario(
        payload: dict[str, Any],
        trusted,
        *,
        comparison: bool = False,
) -> dict[str, Any]:
    try:
        request = ScenarioRequest.from_dict(payload)
    except ValueError as exc:
        raise HTTPException(status_code=422, detail=str(exc)) from exc

    try:
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
        runner = execute_comparison if comparison else execute
        result = await run_in_threadpool(runner, request, isolated_trusted)
        await run_in_threadpool(repository.complete, result)
        return result
    except (OSError, RuntimeError, ValueError) as exc:
        logger.exception("Scenario execution failed")
        try:
            await run_in_threadpool(repository.release, request.scenario_id)
        except Exception:
            logger.exception("Could not release failed scenario reservation")
        raise HTTPException(status_code=502, detail="Scenario execution failed") from exc



@app.post("/api/v1/networks/{network_id}/comparisons")
async def create_network_comparison(
        network_id: str,
        payload: dict[str, Any],
) -> dict[str, Any]:
    try:
        registry = load_network_registry()
    except (RuntimeError, OSError, ValueError) as exc:
        logger.exception("Network registry unavailable")
        raise HTTPException(
            status_code=503,
            detail="Network registry unavailable",
        ) from exc

    network = registry.get(network_id)
    if network is None:
        raise HTTPException(status_code=404, detail="Network not found")

    try:
        request = ScenarioRequest.from_dict(payload)

        if request.scenario_type != "capacity_planning":
            raise ValueError("Comparison requires capacity_planning")

        if request.include_layer1:
            raise ValueError("Fast comparison does not support Layer-1")

        if request.demand_model not in {
            "uniform", "bimodal", "gaussian"
        }:
            raise ValueError("Unsupported fast-comparison demand model")

    except ValueError as exc:
        raise HTTPException(status_code=422, detail=str(exc)) from exc

    try:
        trusted = await run_in_threadpool(
            trusted_run_for_network, network
        )
    except (OSError, RuntimeError, ValueError) as exc:
        logger.exception("Network artifact verification failed")
        raise HTTPException(
            status_code=503,
            detail="Network unavailable",
        ) from exc

    return await _run_scenario(
        payload,
        trusted,
        comparison=True,
    )


@app.get("/api/v1/networks")
def list_networks() -> dict[str, Any]:
    try:
        registry = load_network_registry()
    except (RuntimeError, OSError, ValueError) as exc:
        logger.exception("Network registry unavailable")
        raise HTTPException(status_code=503, detail="Network registry unavailable") from exc
    return {"networks": [
        {"id": network.network_id, "name": network.display_name, "solver": network.solver_label}
        for network in registry.values()
    ]}


@app.post("/api/v1/networks/{network_id}/scenarios")
async def create_network_scenario(network_id: str, payload: dict[str, Any]) -> dict[str, Any]:
    try:
        registry = load_network_registry()
    except (RuntimeError, OSError, ValueError) as exc:
        logger.exception("Network registry unavailable")
        raise HTTPException(status_code=503, detail="Network registry unavailable") from exc
    network = registry.get(network_id)
    if network is None:
        raise HTTPException(status_code=404, detail="Network not found")
    try:
        trusted = await run_in_threadpool(trusted_run_for_network, network)
    except (OSError, RuntimeError, ValueError) as exc:
        logger.exception("Network artifact verification failed")
        raise HTTPException(status_code=503, detail="Network unavailable") from exc
    return await _run_scenario(payload, trusted)


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


@app.post("/api/v1/scenarios/{scenario_id}/ask")
async def ask_scenario(scenario_id: str, payload: ScenarioQuestion) -> dict[str, Any]:
    repository = ScenarioRepository(
        Path(os.environ.get("ER_RESPONSE_DIR", "results/scenario_responses"))
    )

    try:
        result = await run_in_threadpool(repository.get, scenario_id)
    except ValueError as exc:
        raise HTTPException(status_code=422, detail=str(exc)) from exc
    except OSError as exc:
        logger.exception("Scenario evidence could not be loaded")
        raise HTTPException(status_code=503, detail="Scenario storage unavailable") from exc

    if result is None:
        raise HTTPException(status_code=404, detail="Scenario not found")
    if result.get("status") != "completed":
        raise HTTPException(status_code=409, detail="Scenario is not completed")
    if result.get("scenario_id") != scenario_id:
        raise HTTPException(status_code=409, detail="Scenario evidence ID mismatch")

    if not os.environ.get("OPENAI_API_KEY"):
        raise HTTPException(status_code=503, detail="Jeylan is not configured")

    from openai import OpenAI

    try:
        return await run_in_threadpool(
            answer_scenario_question,
            result,
            payload.question,
            OpenAI(timeout=30.0, max_retries=0),
            os.environ.get("JEYLAN_MODEL", "gpt-5.4-mini"),
        )
    except ValueError as exc:
        raise HTTPException(status_code=422, detail=str(exc)) from exc
    except Exception as exc:
        logger.exception("Jeylan scenario interpretation failed")
        raise HTTPException(status_code=502, detail="Jeylan interpretation unavailable") from exc
