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
    """Validate a C++ scenario report and derive display data.

    Supports capacity upgrades and single-link failures.
    Does not perform any new routing computations.
    """
    if not isinstance(report, dict):
        raise ValueError("Scenario report must be an object")

    version = report.get("schema_version")
    if version not in {"1.0", "1.1"}:
        raise ValueError("Unsupported scenario report schema version")

    kind = report.get("scenario_kind")
    if kind not in {"capacity_upgrade", "link_failure"}:
        raise ValueError("Unsupported scenario kind")

    solvers = report.get("solver_results")
    if not isinstance(solvers, list) or not solvers:
        raise ValueError("Missing solver_results")

    result = {
        "schema_version": version,
        "scenario_kind": kind,
        "execution_mode": report.get("execution_mode"),
        "original_fingerprint": str(
            report["original_fingerprint"]
        ),
        "modified_fingerprint": str(
            report["modified_fingerprint"]
        ),
        "capacity_upgrades": report.get(
            "capacity_upgrades", []
        ),
        "link_failures": report.get(
            "link_failures", []
        ),
        "solver_results": [],
    }

    if (
            not isinstance(result["capacity_upgrades"], list)
            or not isinstance(result["link_failures"], list)
    ):
        raise ValueError("Invalid scenario modification lists")

    def validate_network(state, name):
        if not isinstance(state, dict):
            raise ValueError(f"Missing {name} network state")

        links = state.get("links")
        if not isinstance(links, list) or not links:
            raise ValueError(f"Missing {name}.links")

        for field in (
                "max_congestion",
                "average_congestion",
                "total_demand_volume",
        ):
            if _require_number(
                    state.get(field), f"{name}.{field}"
            ) < 0:
                raise ValueError(f"Negative {name}.{field}")

        edge_ids = set()

        for link in links:
            if not isinstance(link, dict):
                raise ValueError(f"Invalid {name} link")

            edge_id = link.get("edge_id")

            if (
                    type(edge_id) is not int
                    or edge_id < 0
                    or edge_id in edge_ids
            ):
                raise ValueError(
                    f"Invalid or duplicate {name} edge ID"
                )

            edge_ids.add(edge_id)

            source = link.get("source")
            target = link.get("target")

            if (
                    type(source) is not int
                    or type(target) is not int
                    or source < 0
                    or target < 0
                    or source == target
            ):
                raise ValueError(
                    f"Invalid {name} physical-link endpoints"
                )

            for field in (
                    "load", "capacity", "utilization"
            ):
                value = _require_number(
                    link.get(field),
                    f"{name}.{field}",
                )

                if value < 0:
                    raise ValueError(
                        f"Negative {name}.{field}"
                    )

            if link["capacity"] <= 0:
                raise ValueError(
                    f"Nonpositive {name} capacity"
                )

            if not math.isclose(
                    link["utilization"],
                    link["load"] / link["capacity"],
                    rel_tol=1e-7,
                    abs_tol=1e-7,
            ):
                raise ValueError(
                    f"Inconsistent {name} utilization"
                )

        actual_max = max(
            link["utilization"] for link in links
        )

        if not math.isclose(
                state["max_congestion"],
                actual_max,
                rel_tol=1e-7,
                abs_tol=1e-7,
        ):
            raise ValueError(
                f"Inconsistent {name} maximum congestion"
            )

        return state

    for solver in solvers:
        if not isinstance(solver, dict):
            raise ValueError("Invalid solver result")

        names = solver.get("solver_algorithms")
        if (
                not isinstance(names, list)
                or not names
                or not all(
            isinstance(name, str) and name
            for name in names
        )
        ):
            raise ValueError("Invalid solver_algorithms")

        baseline = validate_network(
            solver.get("baseline"), "baseline"
        )

        reoptimized = validate_network(
            solver.get("reoptimized"), "reoptimized"
        )

        fixed = solver.get("fixed_routing")

        if kind == "capacity_upgrade":
            fixed = validate_network(
                fixed, "fixed_routing"
            )
        elif fixed is not None:
            raise ValueError(
                "Link failure must not have fixed routing"
            )

        b = {
            link["edge_id"]: link
            for link in baseline["links"]
        }

        r = {
            link["edge_id"]: link
            for link in reoptimized["links"]
        }

        f = (
            {
                link["edge_id"]: link
                for link in fixed["links"]
            }
            if fixed is not None else None
        )

        if kind == "capacity_upgrade":
            if b.keys() != r.keys() or b.keys() != f.keys():
                raise ValueError(
                    "Capacity upgrade changed physical links"
                )
        else:
            failed = set(result["link_failures"])

            if (
                    not failed
                    or any(type(e) is not int for e in failed)
                    or not failed <= b.keys()
                    or r.keys() != b.keys() - failed
            ):
                raise ValueError(
                    "Link failure topology mismatch"
                )

        if not math.isclose(
                baseline["total_demand_volume"],
                reoptimized["total_demand_volume"],
                rel_tol=1e-8,
                abs_tol=1e-8,
        ):
            raise ValueError(
                "Scenario changed total demand volume"
            )

        if fixed is not None:
            if not math.isclose(
                    baseline["total_demand_volume"],
                    fixed["total_demand_volume"],
                    rel_tol=1e-8,
                    abs_tol=1e-8,
            ):
                raise ValueError(
                    "Fixed routing changed total demand"
                )

        link_deltas = [
            {
                "edge_id": edge_id,
                "baseline_load": b[edge_id]["load"],
                "fixed_load": (
                    f[edge_id]["load"] if f else None
                ),
                "reoptimized_load": r[edge_id]["load"],
                "delta_from_baseline":
                    r[edge_id]["load"] - b[edge_id]["load"],
                "delta_from_fixed": (
                    r[edge_id]["load"] - f[edge_id]["load"]
                    if f else None
                ),
            }
            for edge_id in sorted(b.keys() & r.keys())
        ]

        states = {
            "baseline": baseline,
            "fixed_routing": fixed,
            "reoptimized": reoptimized,
        }

        bottlenecks = {
            key: (
                sorted(
                    state["links"],
                    key=lambda link: link["utilization"],
                    reverse=True,
                )[:5]
                if state is not None else []
            )
            for key, state in states.items()
        }

        result["solver_results"].append({
            "solver_algorithms": names,
            "baseline_solver_seconds": solver.get(
                "baseline_solver_seconds"
            ),
            "modified_solver_seconds": solver.get(
                "modified_solver_seconds"
            ),
            "baseline_evaluation_seconds": solver.get(
                "baseline_evaluation_seconds"
            ),
            "modified_evaluation_seconds": solver.get(
                "modified_evaluation_seconds"
            ),
            "baseline": baseline,
            "fixed_routing": fixed,
            "reoptimized": reoptimized,
            "comparison": solver.get("comparison"),
            "fixed_routing_comparison": solver.get(
                "fixed_routing_comparison"
            ),
            "performance_metrics": solver.get(
                "performance_metrics"
            ),
            "bottlenecks": bottlenecks,
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
