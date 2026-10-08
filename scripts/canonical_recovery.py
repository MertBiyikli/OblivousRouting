
"""Convert E-Routing Layer-2 records into canonical evidence."""

import math


NONNEGATIVE_FIELDS = (
    "total_demand",
    "unroutable_demand",
    "unroutable_demand_fraction",
    "baseline_congestion",
    "post_failure_congestion",
    "congestion_increase_factor",
    "recomputation_runtime_microseconds",
)


def measured(value):
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        return None

    return value if math.isfinite(value) and value >= 0 else None


def normalize_recovery(evaluation, result):
    """Populate N-1 events, keeping legacy-only outputs partial."""

    summary = evaluation.get("failure_recovery")
    raw_events = evaluation.get("failure_recovery_events")
    n1 = result["n1_recovery"]

    if summary is None or raw_events is None:
        return

    if not isinstance(raw_events, list):
        raise ValueError("failure_recovery_events must be an array")

    tested = summary["tested_links"]

    if len(raw_events) != tested:
        raise ValueError(
            f"N-1 event count {len(raw_events)} != tested_links {tested}"
        )

    events = []

    for raw in raw_events:
        event = {
            "failed_edge_id": raw["failed_edge_id"],
            "source": raw["source"],
            "target": raw["target"],
            "graph_disconnected": raw["graph_disconnected"],
            "recomputation_attempted": raw["recomputation_attempted"],
            "recomputation_succeeded": raw["recomputation_succeeded"],
        }

        event.update({
            field: measured(raw.get(field))
            for field in NONNEGATIVE_FIELDS
        })

        if event["graph_disconnected"] and event["recomputation_attempted"]:
            raise ValueError(
                "Disconnected failure cannot attempt recomputation under current N-1 policy"
            )

        if event["recomputation_succeeded"] and not event["recomputation_attempted"]:
            raise ValueError(
                "Successful recomputation must have been attempted"
            )

        events.append(event)

    ids = [e["failed_edge_id"] for e in events]

    if len(ids) != len(set(ids)):
        raise ValueError("Duplicate failed physical-link IDs")

    disconnected = sum(e["graph_disconnected"] for e in events)
    successful = sum(e["recomputation_succeeded"] for e in events)
    failed = sum(
        e["recomputation_attempted"] and not e["recomputation_succeeded"]
        for e in events
    )

    if (disconnected, successful, failed) != (
            summary["disconnected_failures"],
            summary["successful_recomputations"],
            summary["failed_recomputations"],
    ):
        raise ValueError(
            "N-1 event outcomes disagree with aggregate summary"
        )

    n1["events"] = events
    n1["tested_link_count"] = tested
    n1["status"] = "complete"

    n1["summary"] = {
        key: measured(value)
        if isinstance(value, (int, float)) and not isinstance(value, bool)
        else value
        for key, value in summary.items()
    }

    result["disconnecting_failures"] = [
        {
            "failed_edge_id": e["failed_edge_id"],
            "source": e["source"],
            "target": e["target"],
            "unroutable_demand_fraction": e["unroutable_demand_fraction"],
        }
        for e in events
        if e["graph_disconnected"]
    ]

    result["topology"]["physical_link_count"] = tested

    volumes = {e["total_demand"] for e in events}

    if len(volumes) == 1:
        result["demand"]["total_demand"] = volumes.pop()

    # Overall analysis_metadata.status stays partial because
    # baseline per-link utilization and capacity evidence is missing.
