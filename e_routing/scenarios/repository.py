
"""Persistence and reservation for E-Routing scenario responses."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any
import os
import tempfile


class ScenarioRepository:
    def __init__(self, directory: Path):
        self.directory = directory.resolve()
        self.directory.mkdir(parents=True, exist_ok=True)

    def _path(self, scenario_id: str) -> Path:
        # Validate the scenario ID and return the corresponding file path.
        if (not isinstance(scenario_id, str)
                or not scenario_id
                or len(scenario_id) > 80
                or not all(c.isascii() and (c.isalnum() or c in "-_") for c in scenario_id)):
            raise ValueError("Invalid scenario ID")
        return self.directory / f"{scenario_id}.json"

    def reserve(self, scenario_id: str) -> None:
        # Reserve a scenario ID for execution, creating a placeholder file.
        path = self._path(scenario_id)

        with path.open("x", encoding="utf-8") as file:
            json.dump({"scenario_id": scenario_id, "status": "running"}, file)

    def complete(self, result: dict[str, Any]) -> None:
        path = self._path(result["scenario_id"])

        if not path.is_file():
            raise RuntimeError("Scenario reservation does not exist")

        if self.get(result["scenario_id"]).get("status") != "running":
            raise RuntimeError("Scenario reservation is not running")

        # Serialize before touching the filesystem.
        serialized = json.dumps(result, indent=2, allow_nan=False)

        temporary = None

        try:
            with tempfile.NamedTemporaryFile(
                    mode="w",
                    encoding="utf-8",
                    dir=self.directory,
                    prefix=f".{path.stem}-",
                    suffix=".tmp",
                    delete=False,
            ) as file:
                temporary = Path(file.name)
                file.write(serialized)
                file.flush()
                os.fsync(file.fileno())

            temporary.replace(path)

        finally:
            if temporary is not None:
                temporary.unlink(missing_ok=True)

    def release(self, scenario_id: str) -> None:
        # Release an incomplete reservation.
        path = self._path(scenario_id)

        if path.is_file() and self.get(scenario_id).get("status") == "running":
            path.unlink()

    def save(self, result: dict[str, Any]) -> None:
        # Preserve compatibility with direct result persistence.
        path = self._path(result["scenario_id"])
        with path.open("x", encoding="utf-8") as file:
            json.dump(result, file, indent=2, allow_nan=False)

    def get(self, scenario_id: str) -> dict[str, Any] | None:
        path = self._path(scenario_id)

        if not path.is_file():
            return None

        return json.loads(path.read_text(encoding="utf-8"))
