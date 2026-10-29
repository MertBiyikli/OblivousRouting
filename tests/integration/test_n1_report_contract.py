#!/usr/bin/env python3
"""Standalone tests for the new N-1 JSON contract; no engine required."""
import copy
import sys
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))
from e_routing.api.n1_report_api import validate_n1_report


def network(links):
    return {
        "max_congestion": max((l["utilization"] for l in links), default=0.0),
        "average_congestion": sum(l["utilization"] for l in links) / len(links) if links else 0.0,
        "total_demand_volume": 10.0,
        "links": links,
    }


def link(edge_id, source, target, load):
    return {"edge_id": edge_id, "source": source, "target": target,
            "capacity": 10.0, "load": load, "utilization": load / 10}


def valid_report():
    baseline = network([link(0, 0, 1, 5.0), link(2, 1, 2, 8.0)])
    after = network([link(2, 1, 2, 9.0)])
    return {
        "schema_version": "1.0", "analysis_kind": "n1_resilience",
        "graph_fingerprint": "2082483772139137733", "solver": "electrical",
        "baseline": baseline,
        "summary": {"tested_links": 2, "successful": 1, "disconnected": 1, "execution_failed": 0},
        "outcomes": [
            {"edge_id": 0, "source": 0, "target": 1, "status": "success", "reoptimized": after,
             "comparison": {"congestion_before": 0.8, "congestion_after": 0.9,
                            "absolute_reduction": -0.1, "relative_reduction": -0.125}},
            {"edge_id": 2, "source": 1, "target": 2, "status": "disconnected",
             "reoptimized": None, "comparison": None},
        ]
    }


class TestN1ReportContract(unittest.TestCase):
    def test_valid(self):
        self.assertEqual(validate_n1_report(valid_report())["summary"]["tested_links"], 2)

    def test_duplicate_failure(self):
        x = valid_report()
        x["outcomes"][1]["edge_id"] = 0
        with self.assertRaises(ValueError): validate_n1_report(x)

    def test_incorrect_count(self):
        x = valid_report()
        x["summary"]["disconnected"] = 0
        with self.assertRaises(ValueError): validate_n1_report(x)

    def test_disconnected_with_metrics(self):
        x = valid_report()
        x["outcomes"][1]["reoptimized"] = x["outcomes"][0]["reoptimized"]
        with self.assertRaises(ValueError): validate_n1_report(x)

    def test_success_missing_edge_removal(self):
        x = valid_report()
        x["outcomes"][0]["reoptimized"] = copy.deepcopy(x["baseline"])
        with self.assertRaises(ValueError): validate_n1_report(x)

    def test_mismatched_fingerprint(self):
        x = valid_report()
        x["graph_fingerprint"] = 1234
        with self.assertRaises(ValueError): validate_n1_report(x)

    def test_demand_inconsistent(self):
        x = valid_report()
        x["outcomes"][0]["reoptimized"]["total_demand_volume"] = 15
        with self.assertRaises(ValueError): validate_n1_report(x)

    def test_endpoint_mismatch(self):
        x = valid_report()
        x["outcomes"][0]["target"] = 3
        with self.assertRaises(ValueError): validate_n1_report(x)


if __name__ == "__main__":
    unittest.main()
