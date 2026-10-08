import json
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from analysis_api import AnalysisAPI


class AnalysisAPITests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.example = Path(__file__).resolve().parent / "fixtures" / "canonical" / "sample.json"
        cls.schema = Path(__file__).resolve().parents[1] / "schema" / "analysis-result.schema.json"

    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.data = json.loads(self.example.read_text())
        self.path = Path(self.temp.name) / "sample.json"
        self._write()
        self.api = AnalysisAPI(self.temp.name, self.schema)

    def _write(self):
        self.path.write_text(json.dumps(self.data))

    def test_summary_and_list_runs(self):
        run = self.data["analysis_metadata"]["run_id"]
        self.assertEqual(self.api.list_runs(), [run])
        summary = self.api.get_analysis_summary(run)
        self.assertEqual(summary["event_count"], 36)
        self.assertEqual(summary["overall_status"], "partial")
        self.assertIsNone(summary["capacity_compliance_rate"])

    def test_failure_impact(self):
        run = self.api.list_runs()[0]
        actual = self.api.get_failure_impact(run, 14)
        self.assertEqual(actual["impact"]["failed_edge_id"], 14)
        self.assertEqual(actual["evidence"]["failed_edge_id"], 14)
        self.assertIsNone(actual["capacity_compliant"])
        with self.assertRaises(KeyError):
            self.api.get_failure_impact(run, 999999)

    def test_ranking_and_ties(self):
        run = self.api.list_runs()[0]
        result = self.api.rank_critical_links(run, "congestion_increase_factor", 5)
        values = [x["value"] for x in result["links"]]
        self.assertEqual(values, sorted(values, reverse=True))
        self.assertEqual(len(values), 5)
        with self.assertRaises(ValueError):
            self.api.rank_critical_links(run, "unknown")

    def test_rejects_invalid_result(self):
        self.data["n1_recovery"]["events"].pop()
        self._write()
        api = AnalysisAPI(self.temp.name, self.schema)
        with self.assertRaises(ValueError):
            api.get_analysis_summary(api.list_runs()[0])

    def test_rejects_duplicate_run_id(self):
        (Path(self.temp.name) / "duplicate.json").write_text(self.path.read_text())
        with self.assertRaises(ValueError):
            AnalysisAPI(self.temp.name, self.schema)


if __name__ == "__main__":
    unittest.main()
