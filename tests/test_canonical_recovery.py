
import json
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from scripts.canonical_recovery import normalize_recovery


class TestCanonicalRecovery(unittest.TestCase):

    def setUp(self):
        self.result = {
            "n1_recovery": {
                "status": "partial",
                "tested_link_count": None,
                "events": [],
            },
            "disconnecting_failures": [],
            "topology": {},
            "demand": {},
        }

        self.events = [
            {
                "failed_edge_id": 0,
                "source": 0,
                "target": 1,
                "graph_disconnected": False,
                "total_demand": 100.0,
                "unroutable_demand": 0.0,
                "unroutable_demand_fraction": 0.0,
                "baseline_congestion": 2.0,
                "post_failure_congestion": 2.5,
                "congestion_increase_factor": 1.25,
                "recomputation_attempted": True,
                "recomputation_succeeded": True,
                "recomputation_runtime_microseconds": 45.0,
            },
            {
                "failed_edge_id": 2,
                "source": 1,
                "target": 2,
                "graph_disconnected": True,
                "total_demand": 100.0,
                "unroutable_demand": 30.0,
                "unroutable_demand_fraction": 0.3,
                "baseline_congestion": 2.0,
                "post_failure_congestion": -1.0,
                "congestion_increase_factor": -1.0,
                "recomputation_attempted": False,
                "recomputation_succeeded": False,
                "recomputation_runtime_microseconds": -1.0,
            },
        ]

        self.evaluation = {
            "failure_recovery": {
                "tested_links": 2,
                "disconnected_failures": 1,
                "successful_recomputations": 1,
                "failed_recomputations": 0,
            },
            "failure_recovery_events": self.events,
        }

    def test_all_events_preserved(self):
        normalize_recovery(self.evaluation, self.result)

        self.assertEqual(len(self.result["n1_recovery"]["events"]), 2)
        self.assertEqual(self.result["n1_recovery"]["status"], "complete")

    def test_disconnecting_failure(self):
        normalize_recovery(self.evaluation, self.result)

        event = self.result["n1_recovery"]["events"][1]

        self.assertTrue(event["graph_disconnected"])
        self.assertFalse(event["recomputation_attempted"])
        self.assertEqual(event["unroutable_demand_fraction"], 0.3)
        self.assertIsNone(event["post_failure_congestion"])
        self.assertIsNone(event["congestion_increase_factor"])
        self.assertIsNone(event["recomputation_runtime_microseconds"])
        self.assertEqual(len(self.result["disconnecting_failures"]), 1)

    def test_rejects_inconsistent_event_count(self):
        self.evaluation["failure_recovery"]["tested_links"] = 3

        with self.assertRaises(ValueError):
            normalize_recovery(self.evaluation, self.result)

    def test_rejects_duplicate_link_ids(self):
        self.events[1]["failed_edge_id"] = 0

        with self.assertRaises(ValueError):
            normalize_recovery(self.evaluation, self.result)

    def test_legacy_output_remains_partial(self):
        del self.evaluation["failure_recovery_events"]

        normalize_recovery(self.evaluation, self.result)

        self.assertEqual(self.result["n1_recovery"]["status"], "partial")
        self.assertEqual(self.result["n1_recovery"]["events"], [])


if __name__ == "__main__":
    unittest.main()
