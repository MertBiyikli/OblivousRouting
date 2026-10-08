#!/usr/bin/env python3
"""Read-only deterministic queries over validated E-Routing canonical results."""
import argparse
import json
from pathlib import Path

from jsonschema import Draft202012Validator, FormatChecker
from validate_analysis_results import validate_semantics

METRICS = {
    "congestion_increase_factor": "congestion_increase_factor",
    "post_failure_congestion": "post_failure_congestion",
    "unroutable_demand_fraction": "unroutable_demand_fraction",
    "recomputation_runtime_microseconds": "recomputation_runtime_microseconds",
}


class AnalysisAPI:
    def __init__(self, results_dir, schema_path):
        self.results_dir = Path(results_dir)
        self.schema_path = Path(schema_path)
        schema = json.loads(self.schema_path.read_text(encoding="utf-8"))
        Draft202012Validator.check_schema(schema)
        self.validator = Draft202012Validator(schema, format_checker=FormatChecker())
        self._paths = {}
        self._cache = {}
        for path in sorted(self.results_dir.rglob("*.json")):
            # Use the declared run_id rather than filename; reject ambiguity.
            run_id = json.loads(path.read_text(encoding="utf-8"))["analysis_metadata"]["run_id"]
            if run_id in self._paths:
                raise ValueError(f"Duplicate run_id {run_id}: {self._paths[run_id]} and {path}")
            self._paths[run_id] = path
        if not self._paths:
            raise ValueError(f"No canonical JSON files found under {self.results_dir}")

    def list_runs(self):
        return sorted(self._paths)

    def _load(self, run_id):
        if run_id not in self._paths:
            raise KeyError(f"Unknown run_id {run_id!r}; available: {', '.join(self.list_runs())}")
        if run_id not in self._cache:
            data = json.loads(self._paths[run_id].read_text(encoding="utf-8"))
            issues = [f"schema: {e.message}" for e in self.validator.iter_errors(data)]
            if not issues:
                issues.extend(validate_semantics(data))
            if issues:
                raise ValueError(f"Invalid analysis {run_id}: " + "; ".join(issues))
            self._cache[run_id] = data
        return self._cache[run_id]

    @staticmethod
    def _evidence(data, edge_id=None):
        ev = {"run_id": data["analysis_metadata"]["run_id"], "schema_version": data["schema_version"]}
        if edge_id is not None:
            ev["failed_edge_id"] = edge_id
        return ev

    def get_analysis_summary(self, run_id):
        d = self._load(run_id)
        n1 = d["n1_recovery"]
        s = n1.get("summary")
        return {
            "evidence": self._evidence(d),
            "topology_id": d["topology"]["topology_id"],
            "solver_id": d["solver"]["id"],
            "demand_model": d["demand"]["model"],
            "overall_status": d["analysis_metadata"]["status"],
            "n1_status": n1["status"],
            "tested_link_count": n1["tested_link_count"],
            "event_count": len(n1["events"]),
            "disconnected_failures": s["disconnected_failures"] if s else None,
            "successful_recomputations": s["successful_recomputations"] if s else None,
            "baseline_max_congestion": d["baseline"]["metrics"]["max_congestion"],
            "maximum_post_failure_congestion": s.get("maximum_post_failure_congestion") if s else None,
            "maximum_congestion_increase_factor": s.get("maximum_congestion_increase_factor") if s else None,
            "capacity_compliance_rate": n1.get("capacity_compliance_rate"),
            "limitations": ["Recovery is a simulated recomputation, not measured network failover latency.",
                            "Capacity compliance is unavailable unless explicitly computed."],
        }

    def get_failure_impact(self, run_id, failed_edge_id):
        d = self._load(run_id)
        event = next((e for e in d["n1_recovery"]["events"] if e["failed_edge_id"] == failed_edge_id), None)
        if event is None:
            raise KeyError(f"No N-1 event for physical edge {failed_edge_id} in run {run_id}")
        return {"evidence": self._evidence(d, failed_edge_id), "impact": dict(event),
                "capacity_compliant": None, "capacity_compliance_note": "Per-link capacity evidence unavailable"}

    def rank_critical_links(self, run_id, metric="congestion_increase_factor", limit=5):
        if metric not in METRICS:
            raise ValueError(f"Unsupported metric {metric!r}; choose from {', '.join(METRICS)}")
        if not isinstance(limit, int) or isinstance(limit, bool) or limit < 1:
            raise ValueError("limit must be a positive integer")
        d = self._load(run_id)
        events = d["n1_recovery"]["events"]
        available = [e for e in events if e[METRICS[metric]] is not None]
        ranked = sorted(available, key=lambda e: (-e[METRICS[metric]], e["failed_edge_id"]))
        return {"evidence": self._evidence(d), "metric": metric, "direction": "descending",
                "eligible_count": len(available), "excluded_unavailable_count": len(events) - len(available),
                "links": [{"rank": i + 1, "failed_edge_id": e["failed_edge_id"], "source": e["source"],
                           "target": e["target"], "value": e[METRICS[metric]],
                           "graph_disconnected": e["graph_disconnected"],
                           "evidence": self._evidence(d, e["failed_edge_id"])}
                          for i, e in enumerate(ranked[:limit])]}


def main():
    p = argparse.ArgumentParser(description="Read-only E-Routing canonical evidence queries")
    p.add_argument("--results", type=Path, default=Path("results/canonical"))
    p.add_argument("--schema", type=Path, default=Path("schema/analysis-result.schema.json"))
    sub = p.add_subparsers(dest="command", required=True)
    sub.add_parser("list-runs")
    summary = sub.add_parser("summary")
    summary.add_argument("run_id")
    failure = sub.add_parser("failure")
    failure.add_argument("run_id")
    failure.add_argument("edge_id", type=int)
    rank = sub.add_parser("rank")
    rank.add_argument("run_id")
    rank.add_argument("--metric", choices=sorted(METRICS), default="congestion_increase_factor")
    rank.add_argument("--limit", type=int, default=5)
    args = p.parse_args()
    api = AnalysisAPI(args.results, args.schema)
    try:
        if args.command == "list-runs":
            output = api.list_runs()
        elif args.command == "summary":
            output = api.get_analysis_summary(args.run_id)
        elif args.command == "failure":
            output = api.get_failure_impact(args.run_id, args.edge_id)
        else:
            output = api.rank_critical_links(args.run_id, args.metric, args.limit)
    except (KeyError, ValueError) as ex:
        p.error(str(ex))
    print(json.dumps(output, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
