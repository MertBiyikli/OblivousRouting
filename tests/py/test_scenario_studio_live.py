
"""Real Scenario Studio -> FastAPI -> C++ integration."""

import json
import os

import pytest
from playwright.sync_api import sync_playwright

from test_scenario_studio_workspace import (
    BASE_URL,
    EXPECTED,
    configure_workspace,
    snapshot,
)


@pytest.mark.integration
def test_real_n1_workspace_execution():
    with sync_playwright() as playwright:
        browser = playwright.chromium.launch(headless=True)

        try:
            page = browser.new_page()
            page.set_default_timeout(5000)

            page.goto(
                f"{BASE_URL}/ui",
                wait_until="domcontentloaded",
                timeout=10000,
            )

            page.wait_for_function(
                """() => {
                    const network = document.getElementById('network');
            
                    return network &&
                        !network.disabled &&
                        [...network.options].some(
                            option => option.value === 'backbone-1221'
                        );
                }""",
                timeout=10000,
            )

            page.locator("#network").select_option("backbone-1221")

            configure_workspace(page)

            page.locator("#analysis-type").select_option(
                "n1_resilience"
            )
            page.locator("#execution-mode").select_option("fast")

            with page.expect_request(
                    lambda request:
                    request.method == "POST"
                    and request.url.endswith("/analyses/n1"),
                    timeout=30000,
            ) as request_info:
                with page.expect_response(
                        lambda response:
                        response.request.method == "POST"
                        and response.url.endswith("/analyses/n1"),
                        timeout=600000,
                ) as response_info:
                    page.locator("#run").click()

            request = request_info.value.post_data_json
            response = response_info.value

            assert request["workspace"] == EXPECTED
            assert response.status == 201

            body = response.json()

            assert body["status"] == "completed"
            assert body["summary"] == {
                "tested_links": 9,
                "successful": 8,
                "disconnected": 1,
                "execution_failed": 0,
            }

            analysis_id = body["analysis_id"]

            report_response = page.request.get(
                f"{BASE_URL}/api/v1/networks/backbone-1221/"
                f"analyses/n1/{analysis_id}"
            )

            assert report_response.status == 200

            report = report_response.json()

            assert report["analysis_kind"] == "n1_resilience"
            assert report["solver"] == "electrical_sketching"

            baseline = report["baseline"]

            assert baseline["total_demand_volume"] == pytest.approx(
                16458.0
            )

            capacities = {
                link["edge_id"]: link["capacity"]
                for link in baseline["links"]
            }

            assert capacities[4] == pytest.approx(30.0)

            assert len(report["outcomes"]) == 9

            # Running the analysis must not modify the
            # user's active workspace configuration.
            assert snapshot(page) == EXPECTED

        finally:
            browser.close()
