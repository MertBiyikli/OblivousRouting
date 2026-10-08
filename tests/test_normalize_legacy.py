
import json
import unittest

from scripts.normalize_legacy import normalize


class TestLegacyAdapter(unittest.TestCase):

    def test_multiple_demand_models(self):
        source = {
            "schema_version": "1.0",
            "graph": {
                "name": "test_network",
                "nodes": 4,
                "edges": 5,
            },
            "configuration": {
                "seed": 42,
                "threads": 1,
            },
            "solver_results": [
                {
                    "solver": "Electrical Flow",
                    "solver_type": "electrical_sketching",
                    "status": "ok",
                    "runtime": {
                        "solve_microseconds": 2000
                    },
                    "demand_evaluations": [
                        {
                            "demand_model": "uniform",
                            "congestion": 0.7,
                            "runtime_microseconds": 100,
                        },
                        {
                            "demand_model": "gravity",
                            "congestion": 0.9,
                            "runtime_microseconds": 110,
                        },
                    ],
                }
            ],
        }

        raw = json.dumps(source).encode()

        records = list(normalize(source, raw, "test.json"))

        self.assertEqual(len(records), 2)

        uniform = records[0][1]
        gravity = records[1][1]

        self.assertEqual(
            uniform["baseline"]["metrics"]["max_congestion"],
            0.7,
        )

        self.assertEqual(
            gravity["baseline"]["metrics"]["max_congestion"],
            0.9,
        )

        self.assertEqual(
            uniform["n1_recovery"]["status"],
            "not_available",
        )

        self.assertIsNone(
            uniform["topology"]["sha256"]
        )

    def test_rejects_unsupported_schema(self):
        with self.assertRaises(ValueError):
            list(
                normalize(
                    {"schema_version": "2.0"},
                    b"{}",
                    "test.json",
                )
            )


if __name__ == "__main__":
    unittest.main()
