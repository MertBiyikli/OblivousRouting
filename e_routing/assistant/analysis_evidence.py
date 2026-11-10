"""Bounded evidence summaries for E-Routing analyses."""

from __future__ import annotations

import math
from typing import Any


def number(value: Any) -> float | None:
    if type(value) not in (int, float):
        return None
    if not math.isfinite(value):
        return None
    return float(value)


def top_links(
        links: list,
        selected_edge_id: int | None,
        limit: int = 12,
) -> dict:
    valid = [
        link for link in links
        if isinstance(link, dict)
           and type(link.get("edge_id")) is int
           and number(link.get("utilization")) is not None
    ]

    ranked = sorted(
        valid,
        key=lambda link: link["utilization"],
        reverse=True,
    )

    fields = (
        "edge_id", "source", "target", "capacity",
        "load", "utilization", "baseline_utilization",
        "utilization_change",
    )

    def project(link: dict) -> dict:
        return {
            key: value
            for key in fields
            if (value := link.get(key)) is not None
               and (
                       type(value) not in (int, float)
                       or math.isfinite(value)
               )
        }

    selected = next(
        (
            link for link in valid
            if link["edge_id"] == selected_edge_id
        ),
        None,
    )

    return {
        "total_links": len(valid),
        "highest_utilization_links": [
            project(link) for link in ranked[:limit]
        ],
        "selected_link": (
            project(selected) if selected else None
        ),
        "ranking_is_complete": len(valid) <= limit,
    }


def summarize_analysis(
        network_id: str,
        analysis_id: str,
        analysis_type: str,
        report: dict,
        selected_edge_id: int | None = None,
) -> dict:
    evidence = {
        "network_id": network_id,
        "analysis_id": analysis_id,
        "analysis_type": analysis_type,
        "selected_edge_id": selected_edge_id,
        "evidence_source": "persisted_analysis_report",
    }

    if analysis_type == "traffic_growth":
        baseline = report["baseline"]
        modified = report["modified"]

        evidence["summary"] = {
            "routing_mode": "fixed_routing",
            "baseline_max_congestion": baseline["max_congestion"],
            "modified_max_congestion": modified["max_congestion"],
            "baseline_demand_volume": baseline["total_demand_volume"],
            "modified_demand_volume": modified["total_demand_volume"],
            "relative_congestion_change":
                report.get("relative_congestion_change"),
            "relative_demand_change":
                report.get("relative_demand_change"),
            "global_multiplier":
                report.get("requested_global_scale"),
        }

        evidence["links"] = top_links(
            modified["links"], selected_edge_id
        )

    elif analysis_type in ("capacity_upgrade", "link_failure"):
        solver = report["solver_results"][0]
        baseline = solver["baseline"]
        modified = solver.get("reoptimized")
        fixed = solver.get("fixed_routing")

        evidence["summary"] = {
            "scenario_kind": report["scenario_kind"],
            "baseline_max_congestion": baseline["max_congestion"],
            "reoptimized_max_congestion": (
                modified.get("max_congestion")
                if modified else None
            ),
            "fixed_routing_max_congestion": (
                fixed.get("max_congestion")
                if fixed else None
            ),
            "link_failures": report.get("link_failures", []),
            "capacity_upgrades": report.get("capacity_upgrades", []),
            "disconnected": modified is None,
        }

        # A failed link may not exist in the modified topology.
        # Provide its baseline evidence separately.
        baseline_link = next(
            (
                link for link in baseline.get("links", [])
                if link.get("edge_id") == selected_edge_id
            ),
            None,
        )

        evidence["selected_baseline_link"] = baseline_link

        evidence["links"] = top_links(
            modified.get("links", [])
            if modified else [],
            selected_edge_id,
        )

    elif analysis_type == "n1_resilience":
        baseline = report["baseline"]
        outcomes = report["outcomes"]

        successful = [
            outcome for outcome in outcomes
            if outcome["status"] == "success"
        ]

        ranked = sorted(
            successful,
            key=lambda outcome:
            outcome["reoptimized"]["max_congestion"],
            reverse=True,
        )

        selected_failure = next(
            (
                outcome for outcome in outcomes
                if outcome["edge_id"] == selected_edge_id
            ),
            None,
        )

        def failure_summary(outcome):
            modified = outcome.get("reoptimized")
            baseline_congestion = baseline["max_congestion"]

            post_congestion = (
                modified["max_congestion"]
                if modified is not None
                else None
            )

            absolute_change = (
                post_congestion - baseline_congestion
                if post_congestion is not None
                else None
            )

            relative_change = (
                absolute_change / baseline_congestion
                if absolute_change is not None
                   and baseline_congestion > 0
                else None
            )

            return {
                "edge_id": outcome["edge_id"],
                "source": outcome["source"],
                "target": outcome["target"],
                "status": outcome["status"],
                "post_failure_max_congestion": post_congestion,
                "absolute_congestion_change": absolute_change,
                "relative_congestion_change": relative_change,
            }