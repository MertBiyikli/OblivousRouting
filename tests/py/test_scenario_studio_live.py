
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
            assert report["solver"] == "expander_hierarchy"

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



@pytest.mark.integration
def test_real_shared_workspace_cross_analysis_and_reset():
    """Real browser -> FastAPI -> C++ shared-workspace regression."""

    network_id = "backbone-1221"
    base = f"{BASE_URL}/api/v1/networks/{network_id}"

    expected_workspace = {
        "capacity_upgrades": [
            {"edge_id": 4, "factor": 1.5},
            {"edge_id": 10, "factor": 1.5},
        ],
        "demand_scale": 1.2,
        "commodity_overrides": [
            {"source": 0, "target": 5, "factor": 2.0},
        ],
    }

    def run_analysis(page, analysis_type, endpoint):
        page.locator("#analysis-type").select_option(analysis_type)

        with page.expect_response(
                lambda response: (
                        response.request.method == "POST"
                        and response.url.endswith(endpoint)
                ),
                timeout=600000,
        ) as response_info:
            page.locator("#run").click()

        response = response_info.value
        request = response.request.post_data_json

        assert response.status in (200, 201), response.text

        body = response.json()
        assert body["status"] == "completed"

        return request, body

    with sync_playwright() as playwright:
        browser = playwright.chromium.launch(headless=True)

        try:
            page = browser.new_page()
            page.set_default_timeout(10000)

            page.goto(
                f"{BASE_URL}/ui",
                wait_until="domcontentloaded",
                timeout=15000,
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
                timeout=15000,
            )

            page.locator("#network").select_option(network_id)

            # --------------------------------------------------
            # 1. Configure the complete shared workspace
            # --------------------------------------------------

            page.locator("#traffic-scale").fill("1.2")

            for edge_id in (4, 10):
                page.locator("#edge").fill(str(edge_id))
                page.locator("#factor").fill("1.5")
                page.locator("#add").click()

            page.locator("#traffic-source").fill("0")
            page.locator("#traffic-target").fill("5")
            page.locator("#traffic-factor").fill("2.0")
            page.locator("#traffic-add").click()

            assert snapshot(page) == expected_workspace

            # --------------------------------------------------
            # 2. Link Failure with full workspace
            # --------------------------------------------------

            page.locator("#analysis-type").select_option(
                "link_failure"
            )
            page.locator("#failed-edge").fill("14")

            failure_request, failure_body = run_analysis(
                page,
                "link_failure",
                "/scenarios/detailed-analysis",
            )

            assert failure_request["scenario_kind"] == "link_failure"
            assert failure_request["link_failures"] == [14]
            assert failure_request["workspace"] == expected_workspace

            failure_report = failure_body["analysis"]
            solver = failure_report["solver_results"][0]

            baseline = solver["baseline"]
            modified = solver["reoptimized"]

            assert baseline["max_congestion"] == pytest.approx(
                194.2875
            )
            # Expected result for the current Expander Hierarchy solver,
            # seed 42, shared workspace, and failed physical edge 14.
            assert modified["max_congestion"] == pytest.approx(
                226.23
            )

            assert baseline["total_demand_volume"] == pytest.approx(
                15310.8
            )
            assert modified["total_demand_volume"] == pytest.approx(
                15310.8
            )

            baseline_capacities = {
                link["edge_id"]: link["capacity"]
                for link in baseline["links"]
            }

            modified_capacities = {
                link["edge_id"]: link["capacity"]
                for link in modified["links"]
            }

            assert baseline_capacities[4] == pytest.approx(30.0)
            assert baseline_capacities[10] == pytest.approx(37.5)

            assert modified_capacities[4] == pytest.approx(30.0)
            assert modified_capacities[10] == pytest.approx(37.5)

            assert 14 in baseline_capacities
            assert 14 not in modified_capacities
            assert len(modified["links"]) == 8

            page.locator("#status").get_by_text(
                "Completed", exact=True
            ).wait_for()

            assert snapshot(page) == expected_workspace

            # --------------------------------------------------
            # 3. N-1 resilience on the same workspace
            # --------------------------------------------------

            page.locator("#analysis-type").select_option(
                "n1_resilience"
            )
            page.locator("#execution-mode").select_option("fast")

            n1_request, n1_body = run_analysis(
                page,
                "n1_resilience",
                "/analyses/n1",
            )

            assert n1_request["execution_mode"] == "fast"
            assert n1_request["workspace"] == expected_workspace

            assert n1_body["summary"] == {
                "tested_links": 9,
                "successful": 8,
                "disconnected": 1,
                "execution_failed": 0,
            }

            n1_id = n1_body["analysis_id"]

            n1_response = page.request.get(
                f"{base}/analyses/n1/{n1_id}"
            )

            assert n1_response.status == 200

            n1_report = n1_response.json()
            n1_baseline = n1_report["baseline"]

            assert n1_baseline["max_congestion"] == pytest.approx(
                194.2875
            )
            assert n1_baseline[
                       "total_demand_volume"
                   ] == pytest.approx(15310.8)

            assert len(n1_report["outcomes"]) == 9

            assert snapshot(page) == expected_workspace

            # --------------------------------------------------
            # 4. Traffic Growth on the same workspace
            # --------------------------------------------------

            traffic_request, traffic_body = run_analysis(
                page,
                "traffic_growth",
                "/analyses/traffic",
            )

            assert traffic_request["global_scale"] == pytest.approx(
                1.2
            )

            assert traffic_request["commodity_overrides"] == [
                {"source": 0, "target": 5, "factor": 2.0}
            ]

            assert traffic_request["capacity_upgrades"] == [
                {"edge_id": 4, "capacity_factor": 1.5},
                {"edge_id": 10, "capacity_factor": 1.5},
            ]

            assert traffic_body["execution_mode"] == "fixed_routing"

            traffic_id = traffic_body["analysis_id"]

            traffic_response = page.request.get(
                f"{base}/analyses/traffic/{traffic_id}"
            )

            assert traffic_response.status == 200

            traffic_report = traffic_response.json()

            traffic_baseline = traffic_report["baseline"]
            traffic_modified = traffic_report["modified"]

            assert traffic_baseline[
                       "total_demand_volume"
                   ] == pytest.approx(12477.0)

            assert traffic_modified[
                       "total_demand_volume"
                   ] == pytest.approx(15310.8)

            for stage in (traffic_baseline, traffic_modified):
                capacities = {
                    link["edge_id"]: link["capacity"]
                    for link in stage["links"]
                }

                assert capacities[4] == pytest.approx(30.0)
                assert capacities[10] == pytest.approx(37.5)

            assert snapshot(page) == expected_workspace

            # --------------------------------------------------
            # 5. Reset Network must clear modifications
            #    and stale analysis results
            # --------------------------------------------------

            page.locator("#reset-network").click()

            reset_workspace = snapshot(page)

            assert reset_workspace["capacity_upgrades"] == []
            assert reset_workspace["commodity_overrides"] == []
            assert reset_workspace["demand_scale"] == pytest.approx(
                1.0
            )

            assert page.locator(
                "#workspace-summary"
            ).inner_text() == (
                       "Original network · No modifications"
                   )

            assert page.locator(
                "#baseline-congestion"
            ).inner_text() == "—"

            assert page.locator(
                "#scenario-congestion"
            ).inner_text() == "—"

            assert page.locator(
                "#congestion-delta"
            ).inner_text() == "—"

            assert page.locator(
                "#scenario-analysis"
            ).is_hidden()

            assert page.locator(
                "#graph line"
            ).count() == 0

            # --------------------------------------------------
            # 6. Fresh execution must use original network
            # --------------------------------------------------

            clean_request, clean_body = run_analysis(
                page,
                "traffic_growth",
                "/analyses/traffic",
            )

            assert clean_request["global_scale"] == pytest.approx(
                1.0
            )
            assert clean_request["capacity_upgrades"] == []
            assert clean_request["commodity_overrides"] == []

            clean_response = page.request.get(
                f"{base}/analyses/traffic/"
                f"{clean_body['analysis_id']}"
            )

            assert clean_response.status == 200

            clean_report = clean_response.json()

            assert clean_report["baseline"][
                       "max_congestion"
                   ] == pytest.approx(177.2333333333333)

            assert clean_report["modified"][
                       "max_congestion"
                   ] == pytest.approx(177.2333333333333)

            assert clean_report["baseline"][
                       "total_demand_volume"
                   ] == pytest.approx(12477.0)

            assert clean_report["modified"][
                       "total_demand_volume"
                   ] == pytest.approx(12477.0)

            clean_capacities = {
                link["edge_id"]: link["capacity"]
                for link in clean_report["baseline"]["links"]
            }

            assert clean_capacities[4] == pytest.approx(20.0)
            assert clean_capacities[10] == pytest.approx(25.0)

            assert snapshot(page) == reset_workspace

        finally:
            browser.close()
