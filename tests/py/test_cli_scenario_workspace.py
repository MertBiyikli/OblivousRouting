
import json
import os
import subprocess
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[2]

ENGINE = Path(
    os.environ.get(
        "ER_CLI_TEST_ENGINE",
        str(ROOT / "build/mac-debug-tests-no-openmp/oblivious_routing"),
    )
).resolve()

GRAPH = ROOT / "experiments/datasets/small/Backbone/1221.lgf"


def run_scenario(output: Path, *, upgrade: bool):
    command = [
        str(ENGINE),
        "analyze-scenario",
        "--analysis", "link-failure",
        "--solver", "electrical",
        "--graph", str(GRAPH),
        "--demand", "uniform",
        "--seed", "42",
        "--demand-scale", "1.3",
        "--commodity-scale", "0:2:1.5",
        "--link-failure", "14",
        "--output", str(output),
    ]

    if upgrade:
        command += ["--capacity-upgrade", "4:1.5"]

    result = subprocess.run(
        command,
        cwd=ROOT,
        capture_output=True,
        text=True,
        timeout=60,
        check=False,
    )

    assert result.returncode == 0, result.stderr
    assert output.is_file()

    return json.loads(output.read_text())


def capacity_on_edge(report, edge_id):
    links = report["solver_results"][0]["baseline"]["links"]

    matches = [
        link["capacity"]
        for link in links
        if link["edge_id"] == edge_id
    ]

    assert len(matches) == 1
    return matches[0]


def test_combined_workspace_capacity_and_demand(tmp_path):
    baseline = run_scenario(
        tmp_path / "without_upgrade.json",
        upgrade=False,
        )

    upgraded = run_scenario(
        tmp_path / "with_upgrade.json",
        upgrade=True,
        )

    # The workspace capacity upgrade must be real,
    # not merely recorded as metadata.
    original_capacity = capacity_on_edge(baseline, 4)
    effective_capacity = capacity_on_edge(upgraded, 4)

    assert effective_capacity == pytest.approx(
        original_capacity * 1.5
    )

    # Capacity modification must not change the demands
    # used by the uniform demand model.
    baseline_volume = baseline["workspace"]["effective_demand_volume"]
    upgraded_volume = upgraded["workspace"]["effective_demand_volume"]

    assert baseline_volume == pytest.approx(upgraded_volume)

    # The failure must preserve effective demand volume.
    for report in (baseline, upgraded):
        solver = report["solver_results"][0]

        assert solver["baseline"]["total_demand_volume"] == pytest.approx(
            solver["reoptimized"]["total_demand_volume"]
        )

        assert report["workspace"]["effective_graph_fingerprint"] == (
            report["original_fingerprint"]
        )

        assert report["modified_fingerprint"] != (
            report["original_fingerprint"]
        )

        modified_edges = {
            link["edge_id"]
            for link in solver["reoptimized"]["links"]
        }

        assert 14 not in modified_edges

    # The report must disclose the requested modification.
    assert upgraded["workspace"]["capacity_upgrades"] == [
        {"edge_id": 4, "factor": 1.5}
    ]

    assert upgraded["workspace"]["commodity_overrides"] == [
        {"source": 0, "target": 2, "factor": 1.5}
    ]

    assert upgraded["workspace"]["demand_scale"] == 1.3
    assert upgraded["schema_version"] == "1.2"
    assert upgraded["scenario_kind"] == "link_failure"



def test_capacity_analysis_three_way_comparison(tmp_path):
    output = tmp_path / "capacity_comparison.json"

    command = [
        str(ENGINE),
        "analyze-scenario",
        "--analysis", "capacity",
        "--solver", "electrical",
        "--graph", str(GRAPH),
        "--demand", "uniform",
        "--seed", "42",
        "--demand-scale", "1.3",
        "--capacity-upgrade", "4:1.5",
        "--output", str(output),
    ]

    completed = subprocess.run(
        command,
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=False,
        timeout=60,
    )

    assert completed.returncode == 0, completed.stderr

    report = json.loads(output.read_text())
    solver = report["solver_results"][0]

    assert report["scenario_kind"] == "capacity_upgrade"

    assert report["workspace"]["capacity_upgrades"] == []
    assert len(report["capacity_upgrades"]) == 1

    assert solver["baseline"] is not None
    assert solver["fixed_routing"] is not None
    assert solver["reoptimized"] is not None

    before = capacity_on_edge(report, 4)

    fixed_links = solver["fixed_routing"]["links"]
    after = next(
        link["capacity"]
        for link in fixed_links
        if link["edge_id"] == 4
    )

    assert after == pytest.approx(before * 1.5)

    baseline_volume = solver["baseline"]["total_demand_volume"]

    assert solver["fixed_routing"]["total_demand_volume"] == (
        pytest.approx(baseline_volume)
    )

    assert solver["reoptimized"]["total_demand_volume"] == (
        pytest.approx(baseline_volume)
    )

    assert solver["performance_metrics"] is not None

    assert report["original_fingerprint"] != (
        report["modified_fingerprint"]
    )
