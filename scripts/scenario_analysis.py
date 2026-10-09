from typing import Any


def analyze_link_failure(analysis: dict[str, Any], edge_id: int) -> dict[str, Any]:
    """Retrieve an evaluated N-1 link failure scenario."""

    recovery = analysis.get("n1_recovery", {})
    events = recovery.get("events", [])

    if isinstance(edge_id, bool) or not isinstance(edge_id, int) or edge_id < 0:
        return {
            "status": "invalid_request",
            "message": "edge_id must be a non-negative integer",
        }

    if not isinstance(events, list):
        return {
            "status": "invalid_data",
            "message": "n1_recovery.events must be an array",
        }

    matching = [
        event for event in events
        if isinstance(event, dict) and event.get("failed_edge_id") == edge_id
    ]

    if not matching:
        return {
            "status": "not_available",
            "scenario_type": "single_link_failure",
            "edge_id": edge_id,
            "message": "No evaluated N-1 recovery event exists for this edge.",
        }

    if len(matching) != 1:
        return {
            "status": "invalid_data",
            "message": "Multiple recovery events found for the same edge.",
        }

    event = matching[0]
    baseline = event.get("baseline_congestion")
    post_failure = event.get("post_failure_congestion")

    increase_factor = event.get("congestion_increase_factor")

    return {
        "status": "completed" if event.get("recomputation_succeeded") else "failed",
        "scenario_type": "single_link_failure",
        "source": "stored_n1_recovery",
        "edge_id": edge_id,
        "source_node": event.get("source"),
        "target_node": event.get("target"),
        "baseline_congestion": baseline,
        "post_failure_congestion": post_failure,
        "congestion_increase_factor": increase_factor,
        "congestion_increase_percent": (
            (increase_factor - 1.0) * 100.0
            if isinstance(increase_factor, (int, float))
            else None
        ),
        "graph_disconnected": event.get("graph_disconnected"),
        "recomputation_attempted": event.get("recomputation_attempted"),
        "recomputation_succeeded": event.get("recomputation_succeeded"),
        "unroutable_demand": event.get("unroutable_demand"),
        "unroutable_demand_fraction": event.get("unroutable_demand_fraction"),
        "recomputation_runtime_microseconds": event.get(
            "recomputation_runtime_microseconds"
        ),
        "evidence": {
            "section": "n1_recovery.events",
            "failed_edge_id": edge_id,
        },
    }



def assess_link_maintenance(analysis: dict[str, Any], edge_id: int) -> dict[str, Any]:
    """Assess planned maintenance using recorded N-1 recovery evidence."""

    failure = analyze_link_failure(analysis, edge_id)

    if failure["status"] != "completed":
        return {
            "status": "insufficient_evidence",
            "scenario_type": "planned_link_maintenance",
            "edge_id": edge_id,
            "maintenance_approved": False,
            "failure_analysis": failure,
            "message": "A successful N-1 recovery evaluation is required.",
        }

    findings = []
    concerns = []
    required_checks = []

    baseline = failure["baseline_congestion"]
    post = failure["post_failure_congestion"]
    increase = failure["congestion_increase_percent"]

    findings.append({
        "type": "congestion",
        "baseline": baseline,
        "post_failure": post,
        "increase_percent": round(increase, 3) if increase is not None else None,
    })

    if failure["graph_disconnected"] is True:
        concerns.append("The failure disconnects the network.")

    if failure["unroutable_demand"] is not None and failure["unroutable_demand"] > 0:
        concerns.append("Some traffic demand becomes unroutable.")

    if increase is not None and increase > 0:
        concerns.append(
            f"Modeled congestion increases by {increase:.2f}%."
        )

    required_checks.extend([
        "Verify actual traffic demand against the modeled demand.",
        "Verify post-failure link utilization and capacity constraints.",
        "Confirm the production routing system can apply the recovery.",
        "Assess whether an additional failure during maintenance is tolerable.",
        "Confirm monitoring, rollback procedures and maintenance-window readiness.",
    ])

    return {
        "status": "assessment_complete",
        "scenario_type": "planned_link_maintenance",
        "source": "stored_n1_recovery",
        "edge_id": edge_id,
        "endpoints": [
            failure["source_node"],
            failure["target_node"],
        ],
        "maintenance_approved": False,
        "approval_status": "requires_operational_validation",
        "findings": findings,
        "concerns": concerns,
        "required_checks": required_checks,
        "recovery": {
            "graph_disconnected": failure["graph_disconnected"],
            "recomputation_succeeded": failure["recomputation_succeeded"],
            "unroutable_demand": failure["unroutable_demand"],
            "recomputation_runtime_microseconds":
                failure["recomputation_runtime_microseconds"],
        },
        "evidence": failure["evidence"],
        "limitations": [
            "The assessment uses an existing single-link failure simulation.",
            "No live network telemetry was evaluated.",
            "Solver recomputation time is not actual network failover time.",
            "Approval requires independent operational validation.",
        ],
    }