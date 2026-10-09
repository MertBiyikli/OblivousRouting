"""Contract tests with fake model/API: no key, provider traffic, or engine required."""
import json
import sys
import unittest
from pathlib import Path
from types import SimpleNamespace

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from jeylan import JeylanAssistant, TOOL_DEFINITIONS


class FakeOutput:
    def __init__(self, type, **kw):
        self.type = type
        self.data = {"type": type, **kw}
        for k, v in kw.items():
            setattr(self, k, v)

    def model_dump(self, exclude_none=True):
        return self.data


class FakeResponse:
    status = "completed"

    def __init__(self, output, output_text=""):
        self.output = output
        self.output_text = output_text


class FakeClient:
    def __init__(self, *responses):
        self.queued = list(responses)
        self.requests = []
        self.responses = self

    def create(self, **kw):
        self.requests.append(kw)
        return self.queued.pop(0)


class FakeAPI:
    def __init__(self):
        self.calls = []

    def get_analysis_summary(self, run_id):
        self.calls.append(("summary", run_id))
        return {"maximum_congestion": 2.4, "solver": "electrical_sketching"}

    def get_failure_impact(self, run_id, edge_id):
        self.calls.append(("failure", run_id, edge_id))
        return {"failed_edge_id": edge_id, "congestion_increase_factor": 1.6}

    def rank_critical_links(self, run_id, metric, limit):
        self.calls.append(("rank", run_id, metric, limit))
        return [{"edge_id": 14, "congestion_increase_factor": 1.6}][:limit]


def call(name, args):
    return FakeResponse([FakeOutput("function_call", name=name, arguments=json.dumps(args), call_id="call_1")])


class JeylanTests(unittest.TestCase):
    def setUp(self):
        self.api = FakeAPI()

    def test_summary_tool_and_evidence(self):
        client = FakeClient(call("get_analysis_summary", {}), FakeResponse([FakeOutput("message")], "Max congestion 2.4 [E1]."))
        agent = JeylanAssistant(self.api, client, "run-A")
        self.assertIn("2.4", agent.ask("What is the max congestion?"))
        self.assertEqual(self.api.calls, [("summary", "run-A")])
        self.assertEqual(agent.last_evidence[0].id, "E1")
        self.assertFalse(client.requests[0]["store"])
        self.assertEqual(client.requests[0]["tool_choice"], "required")
        self.assertEqual(client.requests[1]["tool_choice"], "auto")
        results = [x for x in client.requests[1]["input"] if x.get("type") == "function_call_output"]
        self.assertEqual(json.loads(results[0]["output"])["run_id"], "run-A")

    def test_failure_and_ranking(self):
        client = FakeClient(call("get_failure_impact", {"failed_edge_id": 14}), FakeResponse([FakeOutput("message")], "Link 14 has a 1.6x impact [E1]."))
        agent = JeylanAssistant(self.api, client, "run-B")
        agent.ask("What if link 14 fails?")
        self.assertEqual(self.api.calls, [("failure", "run-B", 14)])
        self.assertEqual(agent._dispatch("rank_critical_links", {"metric": "congestion_increase_factor", "limit": 5})[0]["edge_id"], 14)

    def test_run_id_never_model_controlled(self):
        agent = JeylanAssistant(self.api, None, "trusted-run")
        with self.assertRaises(ValueError):
            agent._dispatch("get_analysis_summary", {"run_id": "injected-run"})
        self.assertEqual(self.api.calls, [])

    def test_bad_tools_and_values_are_denied(self):
        agent = JeylanAssistant(self.api, None, "run-A")
        with self.assertRaisesRegex(ValueError, "Unknown"):
            agent._dispatch("execute_shell", {"command": "rm -rf"})
        with self.assertRaises(ValueError):
            agent._dispatch("get_failure_impact", {"failed_edge_id": True})
        with self.assertRaises(ValueError):
            agent._dispatch("rank_critical_links", {"metric": "other", "limit": 5})
        with self.assertRaises(ValueError):
            agent._dispatch("rank_critical_links", {"metric": "congestion_increase_factor", "limit": 1000})
        self.assertEqual(self.api.calls, [])

    def test_error_tool_is_not_evidence(self):
        client = FakeClient(call("execute_shell", {}), FakeResponse([FakeOutput("message")], "Unsupported."))
        agent = JeylanAssistant(self.api, client, "run-A")
        with self.assertRaisesRegex(RuntimeError, "did not retrieve"):
            agent.ask("Run an unsupported command")
        self.assertEqual(agent.history, [])

    def test_multi_turn_and_reset(self):
        client = FakeClient(
            call("get_analysis_summary", {}), FakeResponse([FakeOutput("message")], "First [E1]."),
            call("get_failure_impact", {"failed_edge_id": 14}), FakeResponse([FakeOutput("message")], "Second [E1]."),
        )
        agent = JeylanAssistant(self.api, client, "run-A")
        agent.ask("Summary?")
        agent.ask("And edge 14?")
        self.assertGreater(len(client.requests[2]["input"]), len(client.requests[0]["input"]))
        agent.reset()
        self.assertEqual(agent.history, [])

    def test_all_schemas_strict(self):
        self.assertEqual(len(TOOL_DEFINITIONS), 3)
        for tool in TOOL_DEFINITIONS:
            self.assertTrue(tool["strict"])
            self.assertEqual(tool["parameters"]["required"], list(tool["parameters"]["properties"].keys()))
            self.assertFalse(tool["parameters"]["additionalProperties"])


if __name__ == "__main__":
    unittest.main()
