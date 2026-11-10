
"""Network-scoped fixed-routing Traffic Growth API."""

from __future__ import annotations

from pathlib import Path

from fastapi import APIRouter

# Re-export existing public model and validator names so
# current callers remain compatible.
from e_routing.api.traffic_models import (
    CommodityOverride,
    TrafficRequest,
    validate_traffic_report,
)
from e_routing.api.traffic_service import (
    execute_traffic_analysis,
    load_traffic_analysis,
)


def create_traffic_router(
        reports_root: Path,
        timeout_seconds: int = 120,
) -> APIRouter:
    router = APIRouter()
    root = Path(reports_root).resolve()

    @router.post(
        "/api/v1/networks/{network_id}/analyses/traffic",
        status_code=201,
    )
    async def create_traffic(
            network_id: str,
            request: TrafficRequest,
    ):
        return await execute_traffic_analysis(
            network_id=network_id,
            request=request,
            reports_root=root,
            timeout_seconds=timeout_seconds,
        )

    @router.get(
        "/api/v1/networks/{network_id}/analyses/traffic/{analysis_id}"
    )
    async def get_traffic(
            network_id: str,
            analysis_id: str,
    ):
        return await load_traffic_analysis(
            network_id=network_id,
            analysis_id=analysis_id,
            reports_root=root,
        )

    return router
