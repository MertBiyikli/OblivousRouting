

"""E-Routing CLI compatibility and validation contract.

Run:
    python3 -m pytest tests/py/test_cli_contract.py -vv
"""

from __future__ import annotations

import json
import os
import subprocess
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[2]

DEFAULT_ENGINE = (
        ROOT / "build/mac-debug-tests-no-openmp/oblivious_routing"
)

ENGINE = Path(
    os.environ.get("ER_CLI_TEST_ENGINE", str(DEFAULT_ENGINE))
).resolve()

GRAPH = (
        ROOT / "experiments/datasets/small/Backbone/1221.lgf"
)


@pytest.fixture(scope="module", autouse=True)
def require_engine():
    if not ENGINE.is_file():
        pytest.fail(f"C++ executable not found: {ENGINE}")


def run_cli(*args: str) -> subprocess.CompletedProcess:
    return subprocess.run(
        [str(ENGINE), *args],
        cwd=ROOT,
        capture_output=True,
        text=True,
        timeout=15,
        check=False,
    )


def expect_rejected(
        args: list[str],
        message: str,
):
    result = run_cli(*args)

    assert result.returncode != 0, (
        f"Expected CLI rejection: {args}\n"
        f"stdout:\n{result.stdout}\n"
        f"stderr:\n{result.stderr}"
    )

    assert message.lower() in result.stderr.lower(), (
        f"Expected diagnostic: {message}\n"
        f"Actual stderr:\n{result.stderr}"
    )


def base_args(command: str = "solve") -> list[str]:
    return [
        command,
        "--solver", "electrical",
        "--graph", str(GRAPH),
    ]


# -----------------------------------------------------
# Command discovery
# -----------------------------------------------------

def test_help():
    result = run_cli("--help")

    assert result.returncode == 0
    assert "solve" in result.stdout
    assert "analyze-n1" in result.stdout


def test_solver_listing():
    result = run_cli("--list-solvers")

    assert result.returncode == 0
    assert "electrical" in result.stdout
    assert "raecke_ckr" in result.stdout


def test_unknown_command():
    expect_rejected(
        ["unknown-command"],
        "Unknown command",
    )


# -----------------------------------------------------
# Required arguments
# -----------------------------------------------------

def test_missing_solver():
    expect_rejected(
        ["solve", "--graph", str(GRAPH)],
        "Missing required option: --solver",
    )


def test_missing_graph():
    expect_rejected(
        ["solve", "--solver", "electrical"],
        "Missing required option: --graph",
    )


def test_invalid_solver():
    expect_rejected(
        base_args()[:2]
        + ["nonexistent_solver"]
        + base_args()[3:],
        "Unknown solver",
        )


def test_invalid_threads():
    expect_rejected(
        base_args() + ["--threads", "0"],
        "--threads must be a positive integer",
        )


def test_negative_seed():
    expect_rejected(
        base_args() + ["--seed", "-1"],
        "--seed must be a non-negative integer",
        )


# -----------------------------------------------------
# Demand validation
# -----------------------------------------------------

@pytest.mark.parametrize("scale", [
    "0",
    "-1",
    "nan",
    "inf",
])
def test_invalid_demand_scale(scale):
    expect_rejected(
        base_args() + ["--demand-scale", scale],
        "--demand-scale must be finite and positive",
        )


@pytest.mark.parametrize("commodity", [
    "0:0:1.5",
    "0:2:-1",
    "0:2:nan",
    "0:2:inf",
    "0:2:abc",
    "0:2:1.5junk",
    "0:2",
])
def test_invalid_commodity_overrides(commodity):
    expect_rejected(
        base_args("analyze-n1")
        + [
            "--demand", "uniform",
            "--commodity-scale", commodity,
            "--output", "unused-contract-test.json",
        ],
        "commodity",
        )


def test_duplicate_commodity():
    expect_rejected(
        base_args("analyze-n1")
        + [
            "--demand", "uniform",
            "--commodity-scale",
            "0:2:1.5", "0:2:2.0",
            "--output", "unused-contract-test.json",
        ],
        "Duplicate commodity override",
        )


def test_solve_rejects_unconsumed_commodity():
    expect_rejected(
        base_args()
        + ["--commodity-scale", "0:2:1.5"],
        "--commodity-scale is supported only",
        )


# -----------------------------------------------------
# Capacity and failure validation
# -----------------------------------------------------

@pytest.mark.parametrize("upgrade", [
    "-1:1.5",
    "4:0.5",
    "4:nan",
    "4:inf",
    "4:abc",
    "4:1.5junk",
    "invalid",
])
def test_invalid_capacity_upgrade(upgrade):
    expect_rejected(
        base_args()
        + ["--capacity-upgrade", upgrade],
        "capacity upgrade",
        )



@pytest.mark.parametrize("edge", [
    "-1",
    "14abc",
    "abc",
    "14.5",
])
def test_invalid_failure_edge(edge):
    expect_rejected(
        base_args()
        + ["--link-failure", edge],
        "Invalid link failure edge ID",
        )


def test_duplicate_failed_edges():
    expect_rejected(
        base_args()
        + ["--link-failure", "14", "14"],
        "Duplicate link failure",
        )


# -----------------------------------------------------
# Command compatibility
# -----------------------------------------------------

def test_n1_requires_one_demand_model():
    expect_rejected(
        base_args("analyze-n1")
        + ["--output", "unused-contract-test.json"],
        "exactly one --demand",
        )


def test_n1_rejects_link_failure_flag():
    expect_rejected(
        base_args("analyze-n1")
        + [
            "--demand", "uniform",
            "--link-failure", "14",
            "--output", "unused-contract-test.json",
        ],
        "does not accept",
        )


def test_traffic_requires_routing_table():
    expect_rejected(
        base_args("analyze-traffic")
        + [
            "--demand", "uniform",
            "--output", "unused-contract-test.json",
        ],
        "--load-routing-table",
        )



def test_n1_rejects_research_mode():
    expect_rejected(
        base_args("analyze-n1")
        + [
            "--demand", "uniform",
            "--execution-mode", "research",
            "--output", "unused-contract-test.json",
        ],
        "not implemented yet",
        )


def scenario_args():
    return base_args("analyze-scenario") + [
        "--analysis", "link-failure",
        "--demand", "uniform",
        "--link-failure", "14",
        "--output", "unused-scenario-contract.json",
    ]


def test_scenario_command_is_recognized():
    result = run_cli(*scenario_args())

    assert result.returncode == 0, result.stderr
    assert "Scenario report written to:" in result.stdout


def test_scenario_requires_analysis_kind():
    args = scenario_args()
    index = args.index("--analysis")
    del args[index:index + 2]
    expect_rejected(
        args,
        "requires a valid --analysis",
    )


def test_scenario_rejects_unknown_analysis():
    args = scenario_args()
    args[args.index("--analysis") + 1] = "unknown"

    expect_rejected(
        args,
        "Unsupported analysis kind",
    )


def test_scenario_rejects_multiple_failures():
    args = scenario_args()
    index = args.index("--link-failure")
    args.insert(index + 2, "16")

    expect_rejected(
        args,
        "exactly one --link-failure",
    )


def test_scenario_rejects_missing_failure():
    args = scenario_args()
    index = args.index("--link-failure")
    del args[index:index + 2]

    expect_rejected(
        args,
        "exactly one --link-failure",
    )


def test_analysis_option_rejected_on_solve():
    expect_rejected(
        base_args() + [
            "--analysis", "link-failure",
        ],
        "--analysis requires analyze-scenario",
        )


def test_scenario_accepts_combined_workspace():
    result = run_cli(
        *scenario_args(),
        "--demand-scale", "1.3",
        "--capacity-upgrade", "4:1.5",
        "--commodity-scale", "0:2:1.5",
        "--execution-mode", "standard",
    )

    assert result.returncode == 0, result.stderr
    assert "Scenario report written to:" in result.stdout



def capacity_args():
    return base_args("analyze-scenario") + [
        "--analysis", "capacity",
        "--demand", "uniform",
        "--capacity-upgrade", "4:1.5",
        "--output", "unused-capacity-contract.json",
    ]


def test_capacity_requires_upgrade():
    args = capacity_args()
    index = args.index("--capacity-upgrade")
    del args[index:index + 2]

    expect_rejected(
        args,
        "requires 1 to 100 capacity upgrades"
    )


def test_capacity_rejects_link_failure():
    expect_rejected(
        capacity_args() + ["--link-failure", "14"],
        "does not accept --link-failure",
        )

def test_capacity_accepts_multiple_upgrades(tmp_path):
    output = tmp_path / "multiple_capacity_upgrades.json"

    args = capacity_args()

    output_index = args.index("--output")
    args[output_index + 1] = str(output)

    args.extend([
        "--capacity-upgrade", "6:2.0",
    ])

    result = run_cli(*args)

    assert result.returncode == 0, (
        f"Capacity analysis failed:\n{result.stderr}"
    )

    assert output.is_file()

    report = json.loads(output.read_text(encoding="utf-8"))

    assert report["scenario_kind"] == "capacity_upgrade"

    upgrades = report["capacity_upgrades"]

    assert len(upgrades) == 2

    assert {item["edge_id"] for item in upgrades} == {
        4, 6
    }

    assert report["solver_results"][0]["baseline"]
    assert report["solver_results"][0]["reoptimized"]



def test_traffic_rejects_missing_routing_table(tmp_path):
    missing_table = tmp_path / "missing_routing_table.bin"
    output = tmp_path / "traffic_report.json"

    assert not missing_table.exists()

    expect_rejected(
        base_args("analyze-traffic")
        + [
            "--demand", "uniform",
            "--load-routing-table", str(missing_table),
            "--capacity-upgrade", "4:1.5",
            "--output", str(output),
        ],
        "Cannot open routing table",
        )

    assert not output.exists()