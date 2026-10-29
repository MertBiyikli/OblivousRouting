"""Opt-in C++ scenario recomputation; mount alongside existing fast endpoints.
Requires the C++ `ER_SCENARIO_REPORT_PATH` hook described in INTEGRATION.md.
"""
from __future__ import annotations
import json
import os
from pathlib import Path
import re
import subprocess
from fastapi import APIRouter, HTTPException
from pydantic import BaseModel, Field
from e_routing.api.config import load_network_registry
from e_routing.api.scenario_report_api import build_analysis_payload

_ID = re.compile(r"^[A-Za-z0-9_-]{1,80}$")

class Upgrade(BaseModel):
    edge_id: int = Field(ge=0)
    capacity_factor: float = Field(ge=1.0, le=100.0, allow_inf_nan=False)

class DetailedRequest(BaseModel):
    scenario_id: str
    demand_model: str = "uniform"
    demand_scale: float = Field(default=1.0, gt=0.0, le=100.0, allow_inf_nan=False)
    capacity_upgrades: list[Upgrade] = Field(default_factory=list)


def detailed_report_path(root: Path, network_id: str, scenario_id: str) -> Path:
    if not _ID.fullmatch(network_id) or not _ID.fullmatch(scenario_id):
        raise HTTPException(422, "Invalid network/scenario ID")
    return root / "detailed" / network_id / scenario_id / "report.json"


def create_detailed_analysis_router(reports_root: Path | str = "results/scenarios") -> APIRouter:
    router = APIRouter()
    root = Path(reports_root).resolve()

    @router.post("/api/v1/networks/{network_id}/scenarios/detailed-analysis")
    def run_detailed_analysis(network_id: str, request: DetailedRequest):
        path = detailed_report_path(root, network_id, request.scenario_id)
        registry = load_network_registry()
        network = registry.get(network_id)
        if network is None:
            raise HTTPException(404, "Unknown registered network")
        network.verify()
        if request.demand_model not in {"uniform", "gravity", "bimodal", "gaussian"}:
            raise HTTPException(422, "Unsupported demand model")
        if not request.capacity_upgrades:
            raise HTTPException(422, "Detailed capacity analysis requires at least one upgrade")
        ids = [u.edge_id for u in request.capacity_upgrades]
        if len(ids) != len(set(ids)):
            raise HTTPException(422, "Duplicate upgrade edge IDs")
        # Server-controlled output paths; never trust a client-supplied path.
        path.parent.mkdir(parents=True, exist_ok=True)
        if path.exists() or (path.parent / "request.json").exists():
            raise HTTPException(409, "Scenario ID already exists")
        executable = os.environ.get("ER_ENGINE")
        if not executable:
            raise HTTPException(503, "ER_ENGINE is not configured")
        engine = Path(executable).resolve()
        if not engine.is_file():
            raise HTTPException(503, "Configured C++ engine not found")
        metadata = {"network_id": network_id, "scenario_id": request.scenario_id,
                    "demand_model": request.demand_model, "demand_scale": request.demand_scale,
                    "capacity_upgrades": [u.model_dump() for u in request.capacity_upgrades]}
        # O_EXCL prevents reusing the same ID across concurrent requests.
        try:
            with (path.parent / "request.json").open("x", encoding="utf-8") as handle:
                json.dump(metadata, handle)
        except FileExistsError:
            raise HTTPException(409, "Scenario ID already exists")
        command = [str(engine), "solve", "--solver", network.solver_label,
                   "--graph", str(network.graph), "--demand", request.demand_model,
                   "--demand-scale", str(request.demand_scale), "--seed", "42",
                   "--capacity-upgrade"]
        command += [f"{u.edge_id}:{u.capacity_factor}" for u in request.capacity_upgrades]
        command += ["--output", str(path.parent / "experiment.json")]
        env = os.environ.copy()
        env["ER_SCENARIO_REPORT_PATH"] = str(path)
        env["ER_SCENARIO_EXECUTION_MODE"] = "standard"
        try:
            completed = subprocess.run(command, env=env, cwd=str(root.parent.parent),
                                       capture_output=True, text=True, timeout=180)
        except subprocess.TimeoutExpired:
            raise HTTPException(504, "Detailed C++ scenario exceeded timeout")
        except OSError as exc:
            raise HTTPException(502, f"Could not launch C++ engine: {exc}") from exc
        if completed.returncode != 0:
            raise HTTPException(502, "C++ detailed scenario failed; inspect server logs")
        if not path.is_file():
            raise HTTPException(502, "C++ did not write detailed report; apply report hook")
        try:
            report = json.loads(path.read_text(encoding="utf-8"))
            payload = build_analysis_payload(report)
        except (ValueError, KeyError, TypeError) as exc:
            raise HTTPException(502, f"Invalid C++ detailed report: {exc}") from exc
        return {"scenario_id": request.scenario_id, "network_id": network_id,
                "status": "completed", "report_available": True,
                "analysis": payload}

    return router
