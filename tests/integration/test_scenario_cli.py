
#!/usr/bin/env python3
"""End-to-end regression tests for the C++ scenario CLI."""

import argparse
import json
import math
import subprocess
from pathlib import Path
import os


import json

def reject_duplicate_keys(pairs):
    result = {}

    for key, value in pairs:
        if key in result:
            raise ValueError(
                f"Duplicate JSON key: {key}"
            )
        result[key] = value

    return result


def load_scenario_report(path):
    with open(path, encoding="utf-8") as f:
        return json.load(
            f,
            object_pairs_hook=reject_duplicate_keys
        )


def close(a, b, *, tol=1e-8):
    return math.isclose(a, b, rel_tol=tol, abs_tol=tol)


def by_edge(network):
    return {link["edge_id"]: link for link in network["links"]}


def run_scenario(binary, graph, report, scenario_args):
    # Ensure we never accidentally validate an old report.
    report.parent.mkdir(parents=True, exist_ok=True)
    report.unlink(missing_ok=True)

    command = [
        str(binary),
        "solve",
        "--solver", "electrical",
        "--graph", str(graph),
        "--demand", "uniform",
        *scenario_args,
    ]

    env = os.environ.copy()
    env["ER_SCENARIO_REPORT_PATH"] = str(report)

    result = subprocess.run(
        command,
        text=True,
        capture_output=True,
        check=False,
        env=env,
    )

    print("Command:", " ".join(command))
    print("Report destination:", report)
    print("Exit code:", result.returncode)
    print("STDOUT:", result.stdout)
    print("STDERR:", result.stderr)

    assert result.returncode == 0, (
        f"Scenario execution failed:\n"
        f"{result.stdout}\n{result.stderr}"
    )

    assert report.is_file(), f"Missing report: {report}"



    return load_scenario_report(report)


def assert_network(network):
    assert network is not None
    assert math.isfinite(network["max_congestion"])
    assert network["max_congestion"] >= 0
    assert math.isfinite(network["average_congestion"])

    links = by_edge(network)
    assert len(links) == len(network["links"])

    for link in links.values():
        assert link["capacity"] > 0
        assert link["load"] >= 0
        assert close(
            link["utilization"],
            link["load"] / link["capacity"],
            )

    assert close(
        network["max_congestion"],
        max(link["utilization"] for link in links.values()),
    )


def test_capacity(binary, graph, report):
    data = run_scenario(
        binary, graph, report,
        ["--capacity-upgrade", "2:1.5"],
    )

    assert data["scenario_kind"] == "capacity_upgrade"
    assert data["link_failures"] == []
    assert len(data["capacity_upgrades"]) == 1

    upgrade = data["capacity_upgrades"][0]
    assert upgrade["edge_id"] == 2
    assert close(upgrade["capacity_factor"], 1.5)
    assert close(upgrade["original_capacity"], 50)
    assert close(upgrade["upgraded_capacity"], 75)

    for solver in data["solver_results"]:
        baseline = solver["baseline"]
        fixed = solver["fixed_routing"]
        modified = solver["reoptimized"]

        for network in (baseline, fixed, modified):
            assert_network(network)

        assert close(
            baseline["total_demand_volume"],
            fixed["total_demand_volume"],
        )
        assert close(
            baseline["total_demand_volume"],
            modified["total_demand_volume"],
        )

        before = by_edge(baseline)
        after = by_edge(fixed)

        assert before.keys() == after.keys()

        for edge_id in before:
            # Fixed routing must preserve the physical loads.
            assert close(
                before[edge_id]["load"],
                after[edge_id]["load"],
            )

            if edge_id != 2:
                assert close(
                    before[edge_id]["capacity"],
                    after[edge_id]["capacity"],
                )

        assert close(before[2]["capacity"], 50)
        assert close(after[2]["capacity"], 75)

        assert close(
            after[2]["utilization"],
            before[2]["utilization"] / 1.5,
            )

        assert (
                fixed["max_congestion"]
                <= baseline["max_congestion"] + 1e-8
        )

        assert solver["fixed_routing_comparison"] is not None

    print("PASS: CapacityUpgrade end-to-end")


def test_failure(binary, graph, report):
    data = run_scenario(
        binary, graph, report,
        ["--link-failure", "14"],
    )

    assert data["scenario_kind"] == "link_failure"
    assert data["link_failures"] == [14]
    assert data["capacity_upgrades"] == []

    for solver in data["solver_results"]:
        baseline = solver["baseline"]
        modified = solver["reoptimized"]

        assert_network(baseline)
        assert_network(modified)

        assert solver["fixed_routing"] is None
        assert solver["fixed_routing_comparison"] is None
        assert solver["performance_metrics"] is None

        before = by_edge(baseline)
        after = by_edge(modified)

        assert 14 in before
        assert 14 not in after

        assert after.keys() == before.keys() - {14}

        for edge_id in after:
            assert close(
                before[edge_id]["capacity"],
                after[edge_id]["capacity"],
            )

        assert close(
            baseline["total_demand_volume"],
            modified["total_demand_volume"],
        )

        assert solver["comparison"] is not None

    print("PASS: LinkFailure end-to-end")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", required=True, type=Path)
    parser.add_argument("--graph", required=True, type=Path)
    parser.add_argument(
        "--report-dir",
        default="results/scenarios",
        type=Path,
    )
    args = parser.parse_args()

    binary = args.binary.resolve()
    graph = args.graph.resolve()
    report_dir = args.report_dir.resolve()

    assert binary.is_file()
    assert graph.is_file()

    test_capacity(
        binary, graph,
        report_dir / "capacity.json",
        )

    test_failure(
        binary, graph,
        report_dir / "failure.json",
        )

    print("PASS: Both scenario workflows")


if __name__ == "__main__":
    main()
