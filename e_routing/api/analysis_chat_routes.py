"""Jeylan questions about persisted network analyses."""

from __future__ import annotations

import json
import logging
import os
import re
from pathlib import Path

from fastapi import APIRouter, HTTPException
from fastapi.concurrency import run_in_threadpool
from pydantic import BaseModel, ConfigDict, Field

from e_routing.api.config import load_network_registry
from e_routing.api.n1_report_api import validate_n1_report
from e_routing.api.traffic_routes import validate_traffic_report
from e_routing.assistant.analysis_evidence import summarize_analysis
from e_routing.assistant.scenario_chat import (
    ANALYSIS_INSTRUCTIONS,
    answer_evidence_question,
)
from e_routing.assistant.unified_analysis_chat import (
    answer_unified_analysis,
)
import asyncio

from fastapi import Request
from e_routing.api.workspace_models import WorkspaceSnapshot
from e_routing.assistant.experiment_chat import (
    answer_capacity_experiment,
)


logger = logging.getLogger(__name__)
router = APIRouter()

_ID = re.compile(r"^[A-Za-z0-9_-]{1,64}$")

_ANALYSIS_TYPES = {
    "capacity_upgrade",
    "link_failure",
    "n1_resilience",
    "traffic_growth",
}


class ExperimentQuestion(BaseModel):
    model_config = ConfigDict(extra="forbid")

    question: str = Field(min_length=1, max_length=2000)

    workspace: WorkspaceSnapshot = Field(
        default_factory=WorkspaceSnapshot
    )

    source_analysis_id: str | None = Field(
        default=None,
        pattern=r"^[A-Za-z0-9_-]{1,64}$",
    )


class AnalysisQuestion(BaseModel):
    model_config = ConfigDict(extra="forbid")

    question: str = Field(min_length=1, max_length=2000)
    selected_edge_id: int | None = Field(default=None, ge=0)


def load_analysis_evidence(
        network_id: str,
        analysis_type: str,
        analysis_id: str,
        selected_edge_id: int | None,
) -> dict:
    if not _ID.fullmatch(network_id):
        raise HTTPException(422, "Invalid network ID")

    if not _ID.fullmatch(analysis_id):
        raise HTTPException(422, "Invalid analysis ID")

    if analysis_type not in _ANALYSIS_TYPES:
        raise HTTPException(422, "Unsupported analysis type")

    registry = load_network_registry()

    if network_id not in registry:
        raise HTTPException(404, "Network not found")

    root = Path("results/scenarios").resolve()

    if analysis_type in ("capacity_upgrade", "link_failure"):
        directory = root / "detailed" / network_id / analysis_id
        report_path = directory / "report.json"
        request_path = directory / "request.json"

    elif analysis_type == "n1_resilience":
        directory = root / "n1" / network_id / analysis_id
        report_path = directory / "n1.json"
        request_path = directory / "request.json"

    else:
        directory = root / "traffic" / network_id / analysis_id
        report_path = directory / "traffic.json"
        request_path = directory / "traffic_request.json"

    if not report_path.is_file() or not request_path.is_file():
        raise HTTPException(404, "Completed analysis not found")

    try:
        metadata = json.loads(request_path.read_text("utf-8"))
        report = json.loads(report_path.read_text("utf-8"))

        if not isinstance(metadata, dict) or not isinstance(report, dict):
            raise ValueError("Invalid stored document")

        if metadata.get("network_id") != network_id:
            raise ValueError("Network identity mismatch")

        saved_id = metadata.get(
            "scenario_id"
            if analysis_type in ("capacity_upgrade", "link_failure")
            else "analysis_id"
        )

        if saved_id != analysis_id:
            raise ValueError("Analysis identity mismatch")

        if analysis_type in ("capacity_upgrade", "link_failure"):
            if metadata.get("scenario_kind") != analysis_type:
                raise ValueError("Analysis kind mismatch")

            if report.get("scenario_kind") != analysis_type:
                raise ValueError("Report kind mismatch")

            if report.get("status", "completed") != "completed":
                raise ValueError("Analysis is not completed")

            # Reuse the existing detailed-report adapter as a
            # structural validation step.
            from e_routing.api.scenario_report_api import (
                build_analysis_payload,
            )

            solvers = report.get("solver_results")

            disconnected = (
                    analysis_type == "link_failure"
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

            if disconnected:
                # The C++ result contains no post-failure routing.
                # Do not invent a modified network.
                if len(solvers) != 1:
                    raise ValueError("Unexpected solver result count")
            else:
                build_analysis_payload(report)

        elif analysis_type == "n1_resilience":
            report = validate_n1_report(report)

        else:
            report = validate_traffic_report(report)

        return summarize_analysis(
            network_id=network_id,
            analysis_id=analysis_id,
            analysis_type=analysis_type,
            report=report,
            selected_edge_id=selected_edge_id,
        )

    except (ValueError, TypeError, KeyError, IndexError, OSError) as exc:
        logger.exception("Stored analysis evidence validation failed")
        raise HTTPException(
            503, "Analysis evidence unavailable"
        ) from exc


@router.post(
    "/api/v1/networks/{network_id}/analyses/"
    "{analysis_type}/{analysis_id}/ask"
)
async def ask_analysis(
        network_id: str,
        analysis_type: str,
        analysis_id: str,
        payload: AnalysisQuestion,
):
    evidence = await run_in_threadpool(
        load_analysis_evidence,
        network_id,
        analysis_type,
        analysis_id,
        payload.selected_edge_id,
    )

    if not os.environ.get("OPENAI_API_KEY"):
        raise HTTPException(503, "Jeylan is not configured")

    from openai import OpenAI

    try:
        return await run_in_threadpool(
            answer_evidence_question,
            evidence,
            payload.question,
            OpenAI(timeout=30.0, max_retries=0),
            os.environ.get("JEYLAN_MODEL", "gpt-5.4-mini"),
            instructions=ANALYSIS_INSTRUCTIONS,
            evidence_source="persisted_analysis_report",
            evidence_id=analysis_id,
        )

    except ValueError as exc:
        raise HTTPException(422, str(exc)) from exc

    except Exception as exc:
        logger.exception("Jeylan analysis interpretation failed")
        raise HTTPException(
            502, "Jeylan interpretation unavailable"
        ) from exc



class UnifiedAnalysisQuestion(BaseModel):
    model_config = ConfigDict(extra="forbid")

    question: str = Field(min_length=1, max_length=2000)

    workspace: WorkspaceSnapshot = Field(
        default_factory=WorkspaceSnapshot
    )


# One autonomous experiment at a time per API worker.
_experiment_lock = asyncio.Lock()


@router.post(
    "/api/v1/networks/{network_id}/assistant/"
    "experiments/capacity/ask"
)
async def ask_capacity_experiment(
        network_id: str,
        payload: ExperimentQuestion,
        request: Request,
):
    # Initial deployment is deliberately local and opt-in.
    # A public service requires real user authentication,
    # authorization, quotas and cross-worker rate limits.
    if (
            os.environ.get("JEYLAN_EXPERIMENTS_ENABLED") != "1"
            or request.client is None
            or request.client.host not in ("127.0.0.1", "::1")
    ):
        raise HTTPException(
            403,
            "Autonomous experiments are disabled",
        )

    if not _ID.fullmatch(network_id):
        raise HTTPException(422, "Invalid network ID")

    registry = await run_in_threadpool(load_network_registry)
    network = registry.get(network_id)

    if network is None:
        raise HTTPException(404, "Unknown registered network")

    await run_in_threadpool(network.verify)

    if not os.environ.get("OPENAI_API_KEY"):
        raise HTTPException(503, "Jeylan is not configured")

    try:
        await asyncio.wait_for(
            _experiment_lock.acquire(),
            timeout=0.1,
        )
    except asyncio.TimeoutError:
        raise HTTPException(
            429,
            "Another Jeylan experiment is running",
        )

    try:
        from openai import OpenAI

        return await run_in_threadpool(
            answer_capacity_experiment,
            question=payload.question,
            network_id=network_id,
            workspace=payload.workspace,
            source_analysis_id=payload.source_analysis_id,
            client=OpenAI(
                timeout=30.0,
                max_retries=0,
            ),
            model=os.environ.get(
                "JEYLAN_MODEL",
                "gpt-5.4-mini",
            ),
        )

    except ValueError as exc:
        raise HTTPException(422, str(exc)) from exc

    except HTTPException:
        raise

    except Exception as exc:
        logger.exception(
            "Jeylan capacity experiment failed"
        )
        raise HTTPException(
            502,
            "Jeylan experiment unavailable; check persisted "
            "reports before retrying",
        ) from exc

    finally:
        _experiment_lock.release()


@router.post(
    "/api/v1/networks/{network_id}/assistant/analyses/ask"
)
async def ask_unified_analysis(
        network_id: str,
        payload: UnifiedAnalysisQuestion,
        request: Request,
):
    # Preserve the existing local, opt-in deployment boundary.
    # This is not a substitute for production authentication.
    if (
            os.environ.get("JEYLAN_EXPERIMENTS_ENABLED") != "1"
            or request.client is None
            or request.client.host not in ("127.0.0.1", "::1")
    ):
        raise HTTPException(
            403,
            "Autonomous analyses are disabled",
        )

    if not _ID.fullmatch(network_id):
        raise HTTPException(422, "Invalid network ID")

    registry = await run_in_threadpool(
        load_network_registry
    )
    network = registry.get(network_id)

    if network is None:
        raise HTTPException(
            404,
            "Unknown registered network",
        )

    await run_in_threadpool(network.verify)

    if not os.environ.get("OPENAI_API_KEY"):
        raise HTTPException(
            503,
            "Jeylan is not configured",
        )

    # Share the lock with capacity experiments.
    # Do not allow one N-1 analysis and a capacity experiment
    # to execute concurrently through the assistant.
    try:
        await asyncio.wait_for(
            _experiment_lock.acquire(),
            timeout=0.1,
        )
    except asyncio.TimeoutError:
        raise HTTPException(
            429,
            "Another Jeylan analysis is running",
        )

    try:
        from openai import OpenAI

        return await answer_unified_analysis(
            question=payload.question,
            network_id=network_id,
            workspace=payload.workspace,
            client=OpenAI(
                timeout=30.0,
                max_retries=0,
            ),
            model=os.environ.get(
                "JEYLAN_MODEL",
                "gpt-5.4-mini",
            ),
            reports_root=Path(
                "results/scenarios"
            ).resolve(),
        )

    except ValueError as exc:
        raise HTTPException(
            422,
            str(exc),
        ) from exc

    except HTTPException:
        raise

    except Exception as exc:
        logger.exception(
            "Jeylan unified analysis failed"
        )
        raise HTTPException(
            502,
            "Jeylan analysis unavailable; check persisted "
            "reports before retrying",
        ) from exc

    finally:
        _experiment_lock.release()
