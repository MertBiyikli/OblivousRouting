import json
from pathlib import Path

import pytest
from e_routing.api.scenario_report_api import (
    build_analysis_payload,
    create_scenario_analysis_router,
)

SAMPLE = (
        Path(__file__).resolve().parents[2]
        / "results"
        / "scenarios"
        / "example_scenario_report.json"
)


def test_example_report():
    raw = json.loads(SAMPLE.read_text())
    result = build_analysis_payload(raw)
    assert len(result['solver_results']) == 3
    assert isinstance(result['original_fingerprint'], str)
    electric = next(s for s in result['solver_results'] if s['solver_algorithms'] == ['electrical_sketching'])
    assert abs(electric['performance_metrics']['reoptimized_congestion'] - 149.02) < 1e-8
    assert len(electric['link_deltas']) == 9
    assert len(electric['bottlenecks']['reoptimized']) == 5


def test_rejects_missing_links():
    raw = json.loads(SAMPLE.read_text())
    raw['solver_results'][0]['baseline']['links'].append(raw['solver_results'][0]['baseline']['links'][0])
    with pytest.raises(ValueError, match='duplicate'):
        build_analysis_payload(raw)


def test_route_integration(tmp_path):
    pytest.importorskip('fastapi')
    from fastapi import FastAPI
    from fastapi.testclient import TestClient
    report = tmp_path / 'scenario.json'
    report.write_bytes(SAMPLE.read_bytes())
    app = FastAPI()
    app.include_router(create_scenario_analysis_router(
        reports_root=tmp_path,
        lookup_report=lambda network_id, scenario_id: report if (network_id, scenario_id) == ('backbone-1221', 'upgrade-001') else None,
    ))
    client = TestClient(app)
    response = client.get('/api/v1/networks/backbone-1221/scenarios/upgrade-001/analysis')
    assert response.status_code == 200
    assert len(response.json()['solver_results']) == 3
    assert client.get('/api/v1/networks/backbone-1221/scenarios/unknown/analysis').status_code == 404
