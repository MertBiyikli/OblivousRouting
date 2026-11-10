from __future__ import annotations

import hashlib
import json
import os
import shutil
import subprocess
from pathlib import Path
from uuid import uuid4

from fastapi import APIRouter, File, Form, HTTPException, UploadFile
from fastapi.concurrency import run_in_threadpool


MAX_GRAPH_BYTES = 8 * 1024 * 1024
PREPARATION_TIMEOUT = 120

router = APIRouter()


def imported_networks_root() -> Path:
    return Path(
        os.environ.get(
            "ER_IMPORTED_NETWORKS_DIR",
            "data/imported_networks",
        )
    ).resolve()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()

    with path.open("rb") as stream:
        for chunk in iter(
                lambda: stream.read(1024 * 1024),
                b"",
        ):
            digest.update(chunk)

    return digest.hexdigest()


def prepare_import(
        directory: Path,
        network_id: str,
        display_name: str,
) -> dict:
    graph_path = directory / "graph.lgf"
    routing_dir = directory / "routing"
    report_path = directory / "solve.json"

    routing_dir.mkdir()

    engine_setting = os.environ.get("ER_ENGINE")
    if not engine_setting:
        raise RuntimeError("ER_ENGINE is not configured")

    engine = Path(engine_setting).resolve(strict=True)

    if not engine.is_file() or not os.access(engine, os.X_OK):
        raise RuntimeError("ER_ENGINE is not executable")

    command = [
        str(engine),
        "solve",
        "--solver", "expander",
        "--graph", str(graph_path),
        "--demand", "uniform",
        "--seed", "42",
        "--threads", "1",
        "--save-routing-table", str(routing_dir),
        "--output", str(report_path),
    ]

    result = subprocess.run(
        command,
        cwd=directory,
        capture_output=True,
        text=True,
        timeout=PREPARATION_TIMEOUT,
        check=False,
    )

    if result.returncode != 0:
        raise ValueError("C++ graph preparation failed")

    report = json.loads(report_path.read_text(encoding="utf-8"))

    evaluations = report.get("solver_results", [])
    if (
            len(evaluations) != 1
            or evaluations[0].get("solver_type") != "expander_hierarchy"
            or evaluations[0].get("status") != "ok"
    ):
        raise ValueError("Unexpected routing preparation report")

    graph = report.get("graph", {})
    nodes = graph.get("nodes")
    edges = graph.get("edges")

    if (
            type(nodes) is not int
            or type(edges) is not int
            or nodes < 2
            or edges < 1
    ):
        raise ValueError("Invalid generated graph metadata")

    binaries = list(routing_dir.glob("*.bin"))
    if len(binaries) != 1:
        raise ValueError("Expected exactly one routing table")

    table = binaries[0]

    metadata = {
        "id": network_id,
        "name": display_name,
        "graph": "graph.lgf",
        "routing_table": f"routing/{table.name}",
        "solver": "expander",
        "graph_sha256": sha256_file(graph_path),
        "routing_table_sha256": sha256_file(table),
        "nodes": nodes,
        "edges": edges,
    }

    # Publish the manifest only after all artifacts exist.
    (directory / "metadata.json").write_text(
        json.dumps(metadata, indent=2),
        encoding="utf-8",
    )

    return metadata


@router.post("/api/v1/networks/import", status_code=201)
async def import_network(
        file: UploadFile = File(...),
        name: str = Form(...),
):
    if not name.strip() or len(name.strip()) > 160:
        raise HTTPException(422, "Invalid network name")

    if not file.filename or not file.filename.lower().endswith(".lgf"):
        raise HTTPException(422, "Only LGF files are supported")

    network_id = "uploaded_" + uuid4().hex
    root = imported_networks_root()

    # A staging directory is never exposed as a registered network.
    staging_root = root / ".staging"
    staging = staging_root / network_id
    destination = root / network_id

    try:
        staging.mkdir(parents=True, exist_ok=False)

        graph_path = staging / "graph.lgf"
        received = 0

        with graph_path.open("xb") as output:
            while chunk := await file.read(1024 * 1024):
                received += len(chunk)

                if received > MAX_GRAPH_BYTES:
                    raise HTTPException(
                        413,
                        "Graph upload exceeds the 8 MiB limit",
                    )

                output.write(chunk)

        if received == 0:
            raise HTTPException(422, "Graph file is empty")

        metadata = await run_in_threadpool(
            prepare_import,
            staging,
            network_id,
            name.strip(),
        )

        # Publish the completed network atomically on the same filesystem.
        # Imported registry discovery ignores .staging.
        staging.rename(destination)

        return {
            "network_id": network_id,
            "name": metadata["name"],
            "nodes": metadata["nodes"],
            "edges": metadata["edges"],
            "solver": metadata["solver"],
            "status": "ready",
        }

    except HTTPException:
        raise
    except ValueError as exc:
        raise HTTPException(
            422,
            "Uploaded graph could not be prepared",
        ) from exc
    except subprocess.TimeoutExpired as exc:
        raise HTTPException(
            504,
            "Graph preparation timed out",
        ) from exc
    except (OSError, RuntimeError, json.JSONDecodeError) as exc:
        raise HTTPException(
            503,
            "Network import service unavailable",
        ) from exc
    finally:
        await file.close()
        shutil.rmtree(staging, ignore_errors=True)