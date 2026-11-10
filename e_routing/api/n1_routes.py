
"""Network-scoped N-1 resilience API."""

from __future__ import annotations

from pathlib import Path
from fastapi import APIRouter

from e_routing.api.n1_service import (
    N1Request,
    execute_n1_analysis,
    load_n1_analysis,
)


def create_n1_router(
        reports_root: Path,
        timeout_seconds: float = 600.0,
) -> APIRouter:
    router = APIRouter()
    root = Path(reports_root).resolve()

    @router.post(
        "/api/v1/networks/{network_id}/analyses/n1",
        status_code=201,
    )
    async def create_n1(
            network_id: str,
            request: N1Request,
    ):
        return await execute_n1_analysis(
            network_id=network_id,
            request=request,
            reports_root=root,
            timeout_seconds=timeout_seconds,
        )

    @router.get(
        "/api/v1/networks/{network_id}/analyses/n1/{analysis_id}"
    )
    async def get_n1(
            network_id: str,
            analysis_id: str,
    ):
        return await load_n1_analysis(
            network_id=network_id,
            analysis_id=analysis_id,
            reports_root=root,
        )

    return router
