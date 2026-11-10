
"""Reusable trusted Traffic Growth execution service."""

from __future__ import annotations

import json
import logging
import math
import os
import re
import subprocess
from pathlib import Path

from fastapi import HTTPException
from fastapi.concurrency import run_in_threadpool

from e_routing.api.config import (
    load_network_registry,
    trusted_run_for_network,
)
from e_routing.api.traffic_models import (
    TrafficRequest,
    validate_traffic_report,
)


logger = logging.getLogger(__name__)

_ID = re.compile(r"^[A-Za-z0-9_-]{1,64}$")


def write_traffic_status(
        folder: Path,
        analysis_id: str,
        state: str,
        message: str | None = None,
) -> None:
    payload = {
        "analysis_id": analysis_id,
        "status": state,
    }
    if message is not None:
        payload["message"] = message

    temporary = folder / "status.pending.json"
    destination = folder / "status.json"

    temporary.write_text(
        json.dumps(payload),
        encoding="utf-8",
    )
    os.replace(temporary, destination)


def record_traffic_failure(
        folder: Path,
        analysis_id: str,
        message: str,
) -> None:
    try:
        write_traffic_status(
            folder,
            analysis_id,
            "failed",
            message,
        )
    except OSError:
        logger.exception(
            "Could not persist failed status for %s",
            analysis_id,
        )


def resolve_traffic_network(network_id: str):
    if not _ID.fullmatch(network_id):
        raise HTTPException(422, "Invalid network ID")

    try:
        registry = load_network_registry()
    except (OSError, ValueError, RuntimeError) as exc:
        raise HTTPException(
            503, "Network registry unavailable"
        ) from exc

    network = registry.get(network_id)
    if network is None:
        raise HTTPException(404, "Network not found")

    try:
        trusted = trusted_run_for_network(network)
    except (OSError, ValueError, RuntimeError) as exc:
        raise HTTPException(
            503, "Network artifacts unavailable"
        ) from exc

    return trusted


def traffic_report_paths(
        reports_root: Path | str,
        network_id: str,
        analysis_id: str,
) -> tuple[Path, Path, Path]:
    if not _ID.fullmatch(network_id):
        raise HTTPException(422, "Invalid network ID")

    if not _ID.fullmatch(analysis_id):
        raise HTTPException(422, "Invalid analysis ID")

    folder = (
            Path(reports_root).resolve()
            / network_id
            / analysis_id
    )

    return (
        folder,
        folder / "traffic.json",
        folder / "traffic.pending.json",
    )


async def execute_traffic_analysis(
        network_id: str,
        request: TrafficRequest,
        reports_root: Path | str,
        timeout_seconds: int = 120,
) -> dict:
    """Run one fixed-routing Traffic Growth analysis.

    The network, executable, graph, and routing table come from
    the trusted registry. The caller provides only validated
    analysis parameters and a server-controlled reports root.
    """

    trusted = await run_in_threadpool(
        resolve_traffic_network,
        network_id,
    )

    folder, report_path, pending_path = traffic_report_paths(
        reports_root,
        network_id,
        request.analysis_id,
    )

    folder.mkdir(parents=True, exist_ok=True)

    reservation = folder / "traffic_request.json"

    try:
        fd = os.open(
            reservation,
            os.O_WRONLY | os.O_CREAT | os.O_EXCL,
            0o600,
            )
    except FileExistsError as exc:
        raise HTTPException(
            409, "Analysis ID already exists"
        ) from exc
    except OSError as exc:
        raise HTTPException(
            503, "Analysis storage unavailable"
        ) from exc

    with os.fdopen(fd, "w", encoding="utf-8") as output:
        json.dump(
            {
                "network_id": network_id,
                **request.model_dump(),
            },
            output,
        )

    try:
        write_traffic_status(
            folder,
            request.analysis_id,
            "running",
        )
    except OSError as exc:
        logger.exception(
            "Unable to persist traffic analysis state"
        )
        raise HTTPException(
            503, "Analysis storage unavailable"
        ) from exc

    command = [
        str(trusted.executable),
        "analyze-traffic",
        "--solver", trusted.solver_label,
        "--graph", str(trusted.graph),
        "--load-routing-table", str(trusted.routing_table),
        "--demand", request.demand_model,
        "--seed", str(request.seed),
        "--demand-scale", str(request.global_scale),
    ]

    if request.commodity_overrides:
        command.append("--commodity-scale")
        command.extend(
            f"{item.source}:{item.target}:{item.factor}"
            for item in request.commodity_overrides
        )

    for upgrade in request.capacity_upgrades:
        command.extend([
            "--capacity-upgrade",
            f"{upgrade.edge_id}:{upgrade.capacity_factor}",
        ])

    command.extend([
        "--output", str(pending_path),
    ])

    def execute():
        return subprocess.run(
            command,
            capture_output=True,
            text=True,
            timeout=timeout_seconds,
            check=False,
        )

    try:
        completed = await run_in_threadpool(execute)

    except subprocess.TimeoutExpired as exc:
        logger.exception("Traffic analysis timed out")
        record_traffic_failure(
            folder,
            request.analysis_id,
            "Traffic engine returned an error",
        )
        raise HTTPException(
            504, "Traffic analysis timed out"
        ) from exc

    except OSError as exc:
        logger.exception("Unable to execute traffic engine")
        record_traffic_failure(
            folder,
            request.analysis_id,
            "Traffic engine returned an error",
        )
        raise HTTPException(
            503, "Traffic engine unavailable"
        ) from exc

    if completed.returncode != 0:
        logger.error(
            "Traffic analysis %s failed (exit=%s)\n"
            "Command: %s\n"
            "stdout: %s\n"
            "stderr: %s",
            request.analysis_id,
            completed.returncode,
            command,
            completed.stdout[-4000:],
            completed.stderr[-4000:],
        )
        record_traffic_failure(
            folder,
            request.analysis_id,
            "Traffic engine returned an error",
        )
        raise HTTPException(
            502, "Traffic analysis failed"
        )

    try:
        report = validate_traffic_report(
            json.loads(
                pending_path.read_text(encoding="utf-8")
            )
        )

        if not math.isclose(
                report["requested_global_scale"],
                request.global_scale,
                rel_tol=1e-10,
        ):
            raise ValueError("Global scale mismatch")

        expected = sorted(
            (
                item.source,
                item.target,
                item.factor,
            )
            for item in request.commodity_overrides
        )

        actual = sorted(
            (
                item["source"],
                item["target"],
                item["factor"],
            )
            for item in report["commodity_overrides"]
        )

        if actual != expected:
            raise ValueError("Commodity overrides mismatch")

        # Publish the completed report only after validation.
        pending_path.replace(report_path)

    except (OSError, ValueError) as exc:
        logger.exception(
            "Traffic analysis report validation failed"
        )
        record_traffic_failure(
            folder,
            request.analysis_id,
            "Traffic engine returned an error",
        )
        raise HTTPException(
            502, "Invalid traffic analysis report"
        ) from exc

    try:
        write_traffic_status(
            folder,
            request.analysis_id,
            "completed",
        )
    except OSError as exc:
        logger.exception(
            "Unable to persist completed state"
        )
        raise HTTPException(
            503, "Analysis storage unavailable"
        ) from exc

    return {
        "network_id": network_id,
        "analysis_id": request.analysis_id,
        "status": "completed",
        "execution_mode": "fixed_routing",
        "summary": {
            "baseline_max_congestion":
                report["baseline"]["max_congestion"],
            "modified_max_congestion":
                report["modified"]["max_congestion"],
            "baseline_demand_volume":
                report["baseline"]["total_demand_volume"],
            "modified_demand_volume":
                report["modified"]["total_demand_volume"],
        },
    }


async def load_traffic_analysis(
        network_id: str,
        analysis_id: str,
        reports_root: Path | str,
) -> dict:
    """Retrieve the stored report or current execution status."""

    await run_in_threadpool(
        resolve_traffic_network,
        network_id,
    )

    folder, report_path, _ = traffic_report_paths(
        reports_root,
        network_id,
        analysis_id,
    )

    reservation = folder / "traffic_request.json"
    status_path = folder / "status.json"

    if not reservation.is_file():
        raise HTTPException(
            404, "Traffic analysis not found"
        )

    if not report_path.is_file():
        if status_path.is_file():
            try:
                state = json.loads(
                    status_path.read_text(encoding="utf-8")
                )
            except (OSError, ValueError) as exc:
                raise HTTPException(
                    500, "Analysis status invalid"
                ) from exc

            if state.get("analysis_id") != analysis_id:
                raise HTTPException(
                    500, "Analysis ID mismatch"
                )

            if state.get("status") in {
                "running", "failed"
            }:
                return state

        return {
            "analysis_id": analysis_id,
            "status": "running",
        }

    try:
        return validate_traffic_report(
            json.loads(
                report_path.read_text(encoding="utf-8")
            )
        )
    except (OSError, ValueError) as exc:
        raise HTTPException(
            500, "Stored traffic report invalid"
        ) from exc
