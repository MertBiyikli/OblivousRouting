
"""Trusted, server-side configuration for scenario execution."""

from __future__ import annotations
import hashlib
import json
import os
from pathlib import Path
import re
from dataclasses import dataclass
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




_ID = re.compile(r"^[A-Za-z0-9_-]{1,80}$")
_SHA = re.compile(r"^[a-f0-9]{64}$")


@dataclass(frozen=True)
class NetworkRegistration:
    network_id: str
    display_name: str
    graph: Path
    routing_table: Path
    solver_label: str
    graph_sha256: str
    routing_table_sha256: str

    def verify(self) -> None:
        """Detect changed/swapped graph and routing artifacts before execution."""
        for path, expected in ((self.graph, self.graph_sha256),
                               (self.routing_table, self.routing_table_sha256)):
            digest = hashlib.sha256()
            with path.open("rb") as file:
                for chunk in iter(lambda: file.read(1024 * 1024), b""):
                    digest.update(chunk)
            if digest.hexdigest() != expected:
                raise RuntimeError(f"Registered artifact digest mismatch for {self.network_id}")


def load_network_registry() -> dict[str, NetworkRegistration]:
    """Load administrator-supplied registrations; never accept client paths."""
    config_path = os.environ.get("ER_NETWORK_REGISTRY")
    if not config_path:
        raise RuntimeError("ER_NETWORK_REGISTRY is not configured")
    path = Path(config_path).resolve(strict=True)
    raw = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(raw, dict) or set(raw) != {"networks"} or not isinstance(raw["networks"], list):
        raise ValueError("Invalid network registry document")
    registry: dict[str, NetworkRegistration] = {}
    fields = {"id", "name", "graph", "routing_table", "solver", "graph_sha256", "routing_table_sha256"}
    for item in raw["networks"]:
        if not isinstance(item, dict) or set(item) != fields:
            raise ValueError("Invalid network registration")
        if any(not isinstance(item[k], str) or not item[k] for k in fields):
            raise ValueError("Invalid network registration values")
        network_id = item["id"]
        if not _ID.fullmatch(network_id) or network_id in registry:
            raise ValueError("Invalid or duplicate network ID")
        if not _SHA.fullmatch(item["graph_sha256"]) or not _SHA.fullmatch(item["routing_table_sha256"]):
            raise ValueError("Invalid artifact digest")
        if len(item["name"]) > 160 or not _ID.fullmatch(item["solver"]):
            raise ValueError("Invalid network name or solver")
        # Relative paths are anchored to the administrator-owned registry file.
        graph = (path.parent / item["graph"]).resolve(strict=True)
        routing_table = (path.parent / item["routing_table"]).resolve(strict=True)
        registry[network_id] = NetworkRegistration(
            network_id=network_id,
            display_name=item["name"],
            graph=graph,
            routing_table=routing_table,
            solver_label=item["solver"],
            graph_sha256=item["graph_sha256"],
            routing_table_sha256=item["routing_table_sha256"],
        )
    return registry


def trusted_run_for_network(network: NetworkRegistration) -> TrustedRun:
    # Check artifacts before the subprocess is invoked; C++ still verifies routing validity.
    network.verify()
    engine = os.environ.get("ER_ENGINE")
    if not engine:
        raise RuntimeError("ER_ENGINE is not configured")
    return TrustedRun(
        executable=Path(engine),
        graph=network.graph,
        routing_table=network.routing_table,
        solver_label=network.solver_label,
        workspace=Path(os.environ.get("ER_WORKSPACE", "results/scenarios")),
    )