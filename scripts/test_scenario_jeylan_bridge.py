import unittest
from pathlib import Path
from unittest.mock import patch

from scenario_jeylan_bridge import ScenarioBridge
from scenario_runtime import TrustedRun


class BridgeTests(unittest.TestCase):
    def setUp(self):
        self.bridge = ScenarioBridge(TrustedRun(Path('engine'), Path('graph'), Path('table'), 'electrical', Path('workspace')))

    def test_rejects_untrusted_paths(self):
        with self.assertRaises(ValueError):
            self.bridge.evaluate({'demand_model': 'uniform', 'include_layer1': True, 'graph': '/tmp/other'})

    def test_rejects_unsupported_model(self):
        with self.assertRaises(ValueError):
            self.bridge.evaluate({'demand_model': 'custom', 'include_layer1': True})

    def test_rejects_nonboolean_flag(self):
        with self.assertRaises(ValueError):
            self.bridge.evaluate({'demand_model': 'uniform', 'include_layer1': 'true'})

    @patch('scenario_jeylan_bridge.execute')
    def test_executes_using_trusted_paths_and_returns_evidence(self, run):
        run.return_value = {'status': 'completed', 'execution_mode': 'persisted_routing_scheme', 'metrics': {'congestion': 186.006, 'layer1': {'tested_links': 9}}}
        result = self.bridge.evaluate({'demand_model': 'uniform', 'include_layer1': True})
        self.assertEqual(result['metrics']['congestion'], 186.006)
        request, trusted = run.call_args.args
        self.assertIs(trusted, self.bridge.trusted)
        self.assertEqual(request.demand_model, 'uniform')
        self.assertTrue(request.include_layer1)
        self.assertTrue(request.scenario_id.startswith('jeylan-'))


if __name__ == '__main__':
    unittest.main()

