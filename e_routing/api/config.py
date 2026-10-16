
"""Trusted, server-side configuration for scenario execution."""

from __future__ import annotations

import os
from pathlib import Path

from e_routing.scenarios.runtime import TrustedRun


def load_trusted_config() -> TrustedRun:
    required = ("ER_ENGINE", "ER_GRAPH", "ER_ROUTING_TABLE", "ER_SOLVER")
    missing = [name for name in required if not os.environ.get(name)]

    if missing:
        raise RuntimeError(f"Missing server configuration: {', '.join(missing)}")

    return TrustedRun(
        executable=Path(os.environ["ER_ENGINE"]),
        graph=Path(os.environ["ER_GRAPH"]),
        routing_table=Path(os.environ["ER_ROUTING_TABLE"]),
        solver_label=os.environ["ER_SOLVER"],
        workspace=Path(os.environ.get("ER_WORKSPACE", "results/scenarios")),
    )
