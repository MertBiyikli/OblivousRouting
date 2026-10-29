
import json
import math
import subprocess
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[2]

ENGINE = ROOT / "build/mac-debug-tests-no-openmp/oblivious_routing"
GRAPH = ROOT / "experiments/datasets/small/Backbone/1221.lgf"
ROUTING = ROOT / "results/routing/backbone-1221/0_electrical_sketching.bin"


def analyze(tmp_path, scale, overrides=()):
    output = tmp_path / f"traffic-{scale}-{len(overrides)}.json"

    command = [
        str(ENGINE),
        "analyze-traffic",
        "--solver", "electrical",
        "--graph", str(GRAPH),
        "--load-routing-table", str(ROUTING),
        "--demand", "uniform",
        "--seed", "42",
        "--demand-scale", str(scale),
    ]

    if overrides:
        command.append("--commodity-scale")
        command.extend(
            f"{source}:{target}:{factor}"
            for source, target, factor in overrides
        )

    command.extend(["--output", str(output)])

    completed = subprocess.run(
        command,
        capture_output=True,
        text=True,
        timeout=60,
    )

    assert completed.returncode == 0, completed.stderr
    assert output.is_file()

    return json.loads(output.read_text())


def indexed_links(report, state):
    return {
        link["edge_id"]: link
        for link in report[state]["links"]
    }


@pytest.fixture(scope="module", autouse=True)
def require_artifacts():
    for path in (ENGINE, GRAPH, ROUTING):
        if not path.is_file():
            pytest.skip(f"Missing engine test artifact: {path}")


@pytest.mark.parametrize("scale", [0.5, 1.0, 1.2, 2.0])
def test_uniform_scaling(tmp_path, scale):
    report = analyze(tmp_path, scale)

    baseline = indexed_links(report, "baseline")
    modified = indexed_links(report, "modified")

    assert baseline.keys() == modified.keys()

    assert math.isclose(
        report["modified"]["total_demand_volume"],
        report["baseline"]["total_demand_volume"] * scale,
        rel_tol=1e-9,
        abs_tol=1e-9,
        )

    for edge_id in baseline:
        before = baseline[edge_id]
        after = modified[edge_id]

        assert math.isclose(
            after["capacity"], before["capacity"],
            rel_tol=0,
            abs_tol=0,
        )

        assert math.isclose(
            after["load"], before["load"] * scale,
            rel_tol=1e-8,
            abs_tol=1e-8,
                           ), f"Load mismatch on edge {edge_id}"

        assert math.isclose(
            after["utilization"],
            before["utilization"] * scale,
            rel_tol=1e-8,
            abs_tol=1e-8,
            ), f"Utilization mismatch on edge {edge_id}"

    assert math.isclose(
        report["modified"]["max_congestion"],
        report["baseline"]["max_congestion"] * scale,
        rel_tol=1e-8,
        abs_tol=1e-8,
        )


def test_identity(tmp_path):
    report = analyze(tmp_path, 1.0)

    assert indexed_links(report, "baseline") == indexed_links(
        report, "modified"
    )


def test_single_commodity_override(tmp_path):
    report = analyze(
        tmp_path,
        1.0,
        overrides=[(0, 4, 2.0)],
    )

    baseline = indexed_links(report, "baseline")
    modified = indexed_links(report, "modified")

    changed = []

    for edge_id, before in baseline.items():
        after = modified[edge_id]

        assert after["capacity"] == before["capacity"]

        if not math.isclose(
                after["load"], before["load"],
                rel_tol=1e-9,
                abs_tol=1e-8,
        ):
            changed.append(edge_id)

    # Specific to the current Backbone 1221 electrical routing table:
    # the 0 -> 4 commodity is carried on physical edge 2.
    assert changed == [2]

    assert math.isclose(
        report["baseline"]["max_congestion"],
        report["modified"]["max_congestion"],
        rel_tol=1e-9,
    )


def test_reproducibility(tmp_path):
    a = analyze(tmp_path / "first", 1.2, [(0, 4, 2.0)])
    b = analyze(tmp_path / "second", 1.2, [(0, 4, 2.0)])

    assert a == b
