
import json
import math
import subprocess
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[2]
EXECUTABLE = (
        ROOT / "build/mac-debug-tests-no-openmp/oblivious_routing"
)
GRAPH = (
        ROOT / "experiments/datasets/small/Backbone/1221.lgf"
)


def solve(tmp_path, failure_analysis: bool):
    output = tmp_path / (
        "with_failure.json"
        if failure_analysis
        else "baseline.json"
    )

    command = [
        str(EXECUTABLE),
        "solve",
        "--solver", "electrical",
        "--graph", str(GRAPH),
        "--demand", "uniform",
        "--seed", "42",
        "--threads", "1",
        "--output", str(output),
    ]

    if failure_analysis:
        command.append("--failure-analysis")

    completed = subprocess.run(
        command,
        cwd=ROOT,
        capture_output=True,
        text=True,
        timeout=120,
    )

    assert completed.returncode == 0, (
            completed.stdout + "\n" + completed.stderr
    )

    assert output.is_file()
    return json.loads(output.read_text())


def first_evaluation(data):
    solvers = data["solver_results"]
    assert len(solvers) == 1

    evaluations = solvers[0]["demand_evaluations"]
    assert len(evaluations) == 1

    return evaluations[0]


@pytest.mark.parametrize("failure_analysis", [False, True])
def test_g9_cli_network_data(tmp_path, failure_analysis):
    evaluation = first_evaluation(
        solve(tmp_path, failure_analysis)
    )

    links = evaluation["link_utilizations"]

    # Backbone 1221 currently has nine physical links.
    assert len(links) == 9

    ids = [link["edge_id"] for link in links]
    assert len(ids) == len(set(ids))

    for link in links:
        assert math.isfinite(link["load"])
        assert math.isfinite(link["capacity"])
        assert math.isfinite(link["utilization"])

        assert link["capacity"] > 0.0
        assert link["load"] >= 0.0

        assert math.isclose(
            link["utilization"],
            link["load"] / link["capacity"],
            rel_tol=1e-8,
            abs_tol=1e-8,
            )

    assert math.isclose(
        evaluation["congestion"],
        max(link["utilization"] for link in links),
        rel_tol=1e-8,
        abs_tol=1e-8,
    )


def test_g9_layer1_exposure(tmp_path):
    evaluation = first_evaluation(solve(tmp_path, True))
    failure = evaluation["failure_analysis"]

    assert failure["tested_links"] == 9
    assert failure["most_critical_edge_id"] == 2
    assert failure["runtime_microseconds"] >= 0.0

    assert failure["maximum_lost_traffic_fraction"] == (
        pytest.approx(0.597178809, rel=1e-5)
    )

    assert failure["average_lost_traffic_fraction"] == (
        pytest.approx(0.271926399, rel=1e-5)
    )

    assert (
            failure["maximum_lost_traffic_fraction"]
            >= failure["average_lost_traffic_fraction"]
            >= 0.0
    )

    assert (
            0 <= failure["traffic_carrying_links"] <= 9
    )

    for threshold in (
            "critical_links_10_percent",
            "critical_links_25_percent",
            "critical_links_50_percent",
    ):
        assert 0 <= failure[threshold] <= 9

    assert (
            failure["critical_links_50_percent"]
            <= failure["critical_links_25_percent"]
            <= failure["critical_links_10_percent"]
    )


def test_g9_layer1_disabled(tmp_path):
    evaluation = first_evaluation(solve(tmp_path, False))
    failure = evaluation["failure_analysis"]

    assert failure["tested_links"] == 0
    assert failure["runtime_microseconds"] == -1.0
