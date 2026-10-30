
"""Scenario Studio cross-functional browser regression tests.

Requires a running E-Routing FastAPI server at ER_TEST_BASE_URL.
"""

import os

import pytest
from playwright.sync_api import sync_playwright


BASE_URL = os.getenv(
    "ER_TEST_BASE_URL",
    "http://127.0.0.1:8000",
)



@pytest.fixture
def page():
    with sync_playwright() as playwright:
        browser = playwright.chromium.launch(headless=True)
        context = browser.new_context()
        page = context.new_page()

        page.set_default_timeout(5000)
        page.set_default_navigation_timeout(10000)

        errors = []

        page.on(
            "pageerror",
            lambda error: errors.append(str(error))
        )

        try:
            response = page.goto(
                f"{BASE_URL}/ui",
                wait_until="domcontentloaded",
                timeout=10000,
            )

            assert response is not None
            assert response.status == 200

            page.locator("#network").wait_for(
                state="attached",
                timeout=5000,
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

            assert not errors, (
                f"Scenario Studio JavaScript errors: {errors}"
            )

            page.locator("#network").select_option(
                "backbone-1221"
            )

            yield page

        finally:
            context.close()
            browser.close()


def configure_workspace(page):
    def inspect(label):
        print(f"\n=== {label} ===")
        print("Mode:", page.locator("#analysis-type").input_value())
        print("Scale:", page.locator("#traffic-scale").input_value())
        print("Overrides:", page.locator("#traffic-overrides").inner_text())
        print("Snapshot:", snapshot(page))
        print("Alert:", page.locator("#alert").inner_text())

    page.locator("#analysis-type").select_option("capacity_upgrade")

    page.locator("#edge").fill("4")
    page.locator("#factor").fill("1.5")
    page.locator("#add").click()
    inspect("AFTER CAPACITY UPGRADE")

    page.locator("#analysis-type").select_option("traffic_growth")
    inspect("AFTER MODE SWITCH")

    page.locator("#traffic-scale").fill("1.3")
    inspect("AFTER SETTING SCALE")

    page.locator("#traffic-source").fill("0")
    page.locator("#traffic-target").fill("2")
    page.locator("#traffic-factor").fill("1.5")
    inspect("BEFORE ADDING COMMODITY")

    page.locator("#traffic-add").click()
    inspect("AFTER ADDING COMMODITY")



def snapshot(page):
    """Read the actual workspace state from the running browser."""

    return page.evaluate(
        """() => currentWorkspaceSnapshot()"""
    )


EXPECTED = {
    "capacity_upgrades": [
        {"edge_id": 4, "factor": 1.5},
    ],
    "demand_scale": 1.3,
    "commodity_overrides": [
        {"source": 0, "target": 2, "factor": 1.5},
    ],
}


def test_workspace_composition(page):
    configure_workspace(page)
    print("Analysis:", page.locator("#analysis-type").input_value())
    print("Traffic scale:", page.locator("#traffic-scale").input_value())
    print("Overrides:", page.locator("#traffic-overrides").inner_text())
    print("Workspace:", snapshot(page))
    print("Alert:", page.locator("#alert").inner_text())

    assert snapshot(page) == EXPECTED

    summary = page.locator("#workspace-summary").inner_text()

    assert "1 capacity upgrades" in summary
    assert "1.3" in summary
    assert "1 commodity overrides" in summary


def test_switching_analysis_preserves_workspace(page):
    configure_workspace(page)

    for analysis_type in (
            "n1_resilience",
            "traffic_growth",
            "link_failure",
            "capacity_upgrade",
            "n1_resilience",
    ):
        page.locator("#analysis-type").select_option(analysis_type)

        assert snapshot(page) == EXPECTED, (
            f"Workspace changed when selecting {analysis_type}"
        )

        assert page.locator("#workspace-summary").is_visible()



def test_n1_receives_complete_workspace(page):
    configure_workspace(page)

    page.locator("#analysis-type").select_option("n1_resilience")
    page.locator("#execution-mode").select_option("fast")

    submitted = []

    def capture_request(route):
        submitted.append(route.request.post_data_json)

        route.fulfill(
            status=503,
            content_type="application/json",
            body='{"detail":"Intentional test interception"}',
        )

    page.route(
        "**/api/v1/networks/*/analyses/n1",
        capture_request,
    )

    # Exactly ONE user click.
    page.locator("#run").click()

    # Wait until the UI has processed the intercepted response.
    page.locator("#alert").wait_for(
        state="visible",
        timeout=5000,
    )

    assert len(submitted) == 1, (
        f"Expected one POST, received {len(submitted)}"
    )

    request = submitted[0]

    assert request["workspace"] == EXPECTED
    assert request["execution_mode"] == "fast"
    assert request["demand_model"] == "uniform"
    assert request["seed"] == 42
    assert request["analysis_id"]

    assert snapshot(page) == EXPECTED



def test_unsupported_combination_preserves_state(page):
    configure_workspace(page)

    page.locator("#analysis-type").select_option("traffic_growth")

    submitted = []

    def capture_request(route):
        submitted.append(route.request.post_data_json)
        route.abort()

    page.route(
        "**/api/v1/networks/*/analyses/traffic",
        capture_request,
    )

    page.locator("#run").click()

    assert snapshot(page) == EXPECTED
    assert submitted == []

    assert "capacity upgrades" in (
        page.locator("#alert").inner_text().lower()
    )


def test_reset_network(page):
    configure_workspace(page)

    page.locator("#analysis-type").select_option("n1_resilience")

    page.locator("#reset-network").click()

    assert snapshot(page) == {
        "capacity_upgrades": [],
        "demand_scale": 1.0,
        "commodity_overrides": [],
    }

    assert page.locator("#edge").input_value() == ""
    assert page.locator("#traffic-scale").input_value() == "1.0"

    assert "No modifications" in (
        page.locator("#workspace-summary").inner_text()
    )

    assert page.locator("#graph").locator("*").count() == 0

    assert page.locator("#scenario-analysis").is_hidden()
    assert page.locator("#status").inner_text() == "Ready"


def test_switching_registered_network_resets_workspace(page):
    configure_workspace(page)

    # Dispatch the actual change event even when only one
    # registered network is available in the test environment.
    page.locator("#network").dispatch_event("change")

    assert snapshot(page) == {
        "capacity_upgrades": [],
        "demand_scale": 1.0,
        "commodity_overrides": [],
    }
