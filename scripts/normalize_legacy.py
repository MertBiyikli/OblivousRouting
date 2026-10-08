
#!/usr/bin/env python3

"""
E-Routing: Legacy Benchmark -> Canonical Analysis Adapter.

Converts the existing benchmark JSON (schema_version 1.0)
into one canonical analysis record per solver/demand combination.

Missing evidence is explicitly marked as unavailable.
No routing calculations are performed by this adapter.
"""

import argparse
import hashlib
import json
import math

from datetime import datetime, timezone
from pathlib import Path
from canonical_recovery import  normalize_recovery

SCHEMA_VERSION = "1.0.0"


def optional_number(value):
    if isinstance(value, bool):
        return None

    if not isinstance(value, (int, float)):
        return None

    if not math.isfinite(value) or value < 0:
        return None

    return value


def normalize(source, raw_bytes, source_path):
    if source.get("schema_version") != "1.0":
        raise ValueError("Unsupported benchmark schema version")

    solvers = source.get("solver_results")

    if not isinstance(solvers, list):
        raise ValueError("Missing solver_results array")

    graph = source.get("graph", {})
    config = source.get("configuration", {})

    source_hash = hashlib.sha256(raw_bytes).hexdigest()

    for solver_index, solver in enumerate(solvers):
        evaluations = solver.get("demand_evaluations", [])

        for demand_index, evaluation in enumerate(evaluations):
            model = evaluation["demand_model"]

            run_id = (
                f"{source_hash[:12]}"
                f"-s{solver_index}"
                f"-d{demand_index}"
            )

            layer1 = evaluation.get("failure_analysis", {})
            layer2 = evaluation.get("failure_recovery", {})

            congestion = optional_number(
                evaluation.get("congestion")
            )

            baseline = {
                "status": "partial" if congestion is not None else "unavailable",
                "metrics": {
                    "max_congestion": congestion,
                    "average_utilization": None,
                    "links": [],
                    "overloaded_link_ids": [],
                    "solver_runtime_microseconds": optional_number(
                        solver.get("runtime", {}).get("solve_microseconds")
                    ),
                    "evaluation_runtime_microseconds": optional_number(
                        evaluation.get("runtime_microseconds")
                    ),
                },
            }

            layer1_available = layer1.get("tested_links", 0) > 0
            layer2_available = layer2.get("tested_links", 0) > 0

            result = {
                "schema_version": SCHEMA_VERSION,

                "analysis_metadata": {
                    "run_id": run_id,
                    "scenario_id": f"{graph.get('name', 'unknown')}-{model}",
                    "created_at": datetime.now(timezone.utc).isoformat(),
                    "status": "partial",
                    "enabled_analyses": [
                        "baseline",
                        *(
                            ["layer1_exposure"]
                            if layer1_available else []
                        ),
                        *(
                            ["n1_recovery"]
                            if layer2_available else []
                        ),
                    ],
                },

                "topology": {
                    "topology_id": graph.get("name", "unknown"),
                    "sha256": None,
                    "source": graph.get("path", source_path),
                    "node_count": graph.get("nodes"),
                    "reported_edge_count": graph.get("edges"),
                    "physical_link_count": None,
                },

                "demand": {
                    "demand_id": model,
                    "sha256": None,
                    "source_type": "unknown",
                    "model": model,
                    "scale_factor": None,
                    "total_demand": None,
                    "unit": "unspecified",
                },

                "solver": {
                    "id": solver.get("solver_type", solver.get("solver")),
                    "name": solver.get("solver"),
                    "status": solver.get("status"),
                    "routing_base": solver.get("routing_base"),
                },

                "baseline": baseline,

                "capacity_headroom": {
                    "status": "not_available",
                    "minimum_headroom": None,
                    "violating_link_ids": [],
                },

                "layer1_exposure": {
                    "status": "partial" if layer1_available else "not_available",
                    "tested_link_count": layer1.get("tested_links"),
                    "links": [],
                },

                "n1_recovery": {
                    "status": "partial" if layer2_available else "not_available",
                    "tested_link_count": layer2.get("tested_links"),
                    "events": [],
                    "service_survivability_rate": None,
                    "capacity_compliance_rate": None,
                },

                "critical_links": [],
                "capacity_violations": [],
                "disconnecting_failures": [],
                "maintenance_scenarios": [],
                "growth_scenarios": [],
                "upgrade_candidates": [],
                "recommendations": [],

                "provenance": {
                    "engine_version": None,
                    "engine_commit": None,
                    "solver_id": solver.get("solver_type"),
                    "solver_version": None,
                    "random_seed": config.get("seed"),
                    "thread_count": config.get("threads"),
                    "determinism": "unknown",
                    "configuration": {
                        "legacy_source_sha256": source_hash,
                        "legacy_source_path": str(source_path),
                        "legacy_timestamp": source.get("timestamp"),
                        "legacy_layer1_summary": layer1,
                        "legacy_layer2_summary": layer2,
                    },
                    "tool_calls": [],
                },
            }

            normalize_recovery(evaluation, result)

            yield run_id, result


def main():
    parser = argparse.ArgumentParser(
        description="Normalize E-Routing benchmark results"
    )

    parser.add_argument("input", type=Path)

    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
    )

    args = parser.parse_args()

    raw = args.input.read_bytes()
    source = json.loads(raw)

    records = list(
        normalize(source, raw, str(args.input))
    )

    args.output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    for run_id, result in records:
        output_path = args.output_dir / f"{run_id}.json"

        output_path.write_text(
            json.dumps(
                result,
                indent=2,
                allow_nan=False,
            ) + "\n"
        )

        print(f"Created: {output_path}")

    print(f"Total canonical records: {len(records)}")


if __name__ == "__main__":
    main()
