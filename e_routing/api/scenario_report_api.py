"""Read-only API adapter for the C++ ScenarioResult JSON report.

Mount into an existing FastAPI app; do not replace its POST scenario route.
Use a trusted scenario repository lookup for report paths.
"""
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Callable

from fastapi import APIRouter, HTTPException


ReportLookup = Callable[[str, str], Path | str | None]


def _require_number(value: object, field: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError(f"Invalid numeric field: {field}")
    return float(value)


def build_analysis_payload(report: dict) -> dict:
    """Validate essential schema and derive presentation-only data."""
    if report.get("schema_version") != "1.0":
        raise ValueError("Unsupported scenario report schema version")
    if not isinstance(report.get("solver_results"), list):
        raise ValueError("Missing solver_results array")

    result = {
        "schema_version": "1.0",
        "scenario_kind": report.get("scenario_kind"),
        "execution_mode": report.get("execution_mode"),
        # JS Number cannot precisely represent all uint64 fingerprints.
        "original_fingerprint": str(report["original_fingerprint"]),
        "modified_fingerprint": str(report["modified_fingerprint"]),
        "capacity_upgrades": report.get("capacity_upgrades", []),
        "solver_results": [],
    }
    for solver in report["solver_results"]:
        names = solver.get("solver_algorithms")
        if not isinstance(names, list) or not names or not all(isinstance(n, str) for n in names):
            raise ValueError("Invalid solver_algorithms")
        states = {}
        for key in ("baseline", "fixed_routing", "reoptimized"):
            state = solver.get(key)
            if not isinstance(state, dict) or not isinstance(state.get("links"), list):
                raise ValueError(f"Missing {key} network state")
            _require_number(state.get("max_congestion"), f"{key}.max_congestion")
            ids = set()
            for link in state["links"]:
                edge_id = link.get("edge_id")
                if type(edge_id) is not int or edge_id < 0 or edge_id in ids:
                    raise ValueError(f"Invalid/duplicate {key} physical edge ID")
                ids.add(edge_id)
                for field in ("load", "capacity", "utilization"):
                    if _require_number(link.get(field), f"{key}.{field}") < 0:
                        raise ValueError(f"Negative {key}.{field}")
            states[key] = state
        metrics = solver.get("performance_metrics")
        if not isinstance(metrics, dict):
            raise ValueError("Missing performance_metrics")
        for field in ("baseline_congestion", "fixed_routing_congestion", "reoptimized_congestion"):
            _require_number(metrics.get(field), field)

        b = {link["edge_id"]: link for link in states["baseline"]["links"]}
        f = {link["edge_id"]: link for link in states["fixed_routing"]["links"]}
        r = {link["edge_id"]: link for link in states["reoptimized"]["links"]}
        if b.keys() != f.keys() or b.keys() != r.keys():
            raise ValueError("Physical link sets differ across states")
        link_deltas = [
            {"edge_id": edge_id,
             "baseline_load": b[edge_id]["load"],
             "fixed_load": f[edge_id]["load"],
             "reoptimized_load": r[edge_id]["load"],
             "delta_from_baseline": r[edge_id]["load"] - b[edge_id]["load"],
             "delta_from_fixed": r[edge_id]["load"] - f[edge_id]["load"]}
            for edge_id in sorted(b)
        ]
        result["solver_results"].append({
            "solver_algorithms": names,
            "baseline_solver_seconds": solver.get("baseline_solver_seconds"),
            "modified_solver_seconds": solver.get("modified_solver_seconds"),
            "baseline_evaluation_seconds": solver.get("baseline_evaluation_seconds"),
            "modified_evaluation_seconds": solver.get("modified_evaluation_seconds"),
            "baseline": states["baseline"],
            "fixed_routing": states["fixed_routing"],
            "reoptimized": states["reoptimized"],
            "performance_metrics": metrics,
            "bottlenecks": {key: sorted(state["links"], key=lambda x: x["utilization"], reverse=True)[:5]
                            for key, state in states.items()},
            "link_deltas": link_deltas,
        })
    return result


def create_scenario_analysis_router(*, reports_root: Path, lookup_report: ReportLookup) -> APIRouter:
    """Provide a lookup from the EXISTING scenario repository, not a client-supplied path."""
    router = APIRouter()
    root = Path(reports_root).resolve()

    @router.get("/api/v1/networks/{network_id}/scenarios/{scenario_id}/analysis")
    def get_scenario_analysis(network_id: str, scenario_id: str) -> dict:
        report_path = lookup_report(network_id, scenario_id)
        if report_path is None:
            raise HTTPException(404, "Scenario analysis not found")
        path = Path(report_path).resolve()
        if not path.is_relative_to(root):
            raise HTTPException(500, "Scenario repository points outside report root")
        if not path.is_file():
            raise HTTPException(404, "Scenario report not generated")
        try:
            # Files must have been generated by the trusted C++ pipeline.
            with path.open(encoding="utf-8") as handle:
                report = json.load(handle)
            if not isinstance(report, dict):
                raise ValueError("Report root is not a JSON object")
            return build_analysis_payload(report)
        except (OSError, ValueError, KeyError, TypeError) as exc:
            raise HTTPException(422, f"Invalid scenario report: {exc}") from exc

    return router
