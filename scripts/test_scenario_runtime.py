import unittest
from e_routing.scenarios.runtime import ScenarioRequest

class ScenarioRequestTests(unittest.TestCase):
    def test_baseline(self):
        request = ScenarioRequest.from_dict({'scenario_id':'test-1','scenario_type':'baseline','demand_model':'uniform','include_layer1':True})
        self.assertEqual(request.demand_model, 'uniform')

    def test_reject_unknown_fields(self):
        with self.assertRaises(ValueError):
            ScenarioRequest.from_dict({'scenario_id':'x','scenario_type':'baseline','demand_model':'uniform','include_layer1':True,'command':'rm -rf /'})

    def test_reject_link_removal(self):
        with self.assertRaises(ValueError):
            ScenarioRequest.from_dict({'scenario_id':'x','scenario_type':'link_failure','demand_model':'uniform','include_layer1':True})

    def test_reject_invalid_boolean(self):
        with self.assertRaises(ValueError):
            ScenarioRequest.from_dict({'scenario_id':'x','scenario_type':'baseline','demand_model':'uniform','include_layer1':1})

if __name__ == '__main__':
    unittest.main()
