
import json

from playwright.sync_api import sync_playwright


BASE_URL = "http://127.0.0.1:8000"
NETWORK_ID = "backbone-1221"


def test_jeylan_composite_rendering_and_reset():
    operations = [
        "traffic_growth",
        "n1_resilience",
        "capacity_upgrade",
        "n1_resilience",
    ]

    steps = [
        {
            "step_index": index,
            "operation": operation,
            "analysis_id": f"jeylan-{index + 1:032x}",
            "status": "completed",
            "summary": (
                {"tested_links": 9, "successful": 8}
                if operation == "n1_resilience"
                else {"baseline_max_congestion": 230.4}
            ),
            "input_workspace": {
                "demand_scale": 1.3,
                "capacity_upgrades": [],
                "commodity_overrides": [],
            },
            "output_workspace": {
                "demand_scale": 1.3,
                "capacity_upgrades": [],
                "commodity_overrides": [],
            },
        }
        for index, operation in enumerate(operations)
    ]

    response_body = {
        "network_id": NETWORK_ID,
        "analysis_id": None,
        "analysis_type": "composite",
        "status": "completed",
        "plan": {
            "kind": "composite",
            "demand_model": "uniform",
            "seed": 42,
            "steps": [
                {"operation": operation}
                for operation in operations
            ],
        },
        "completed_steps": steps,
        "comparisons": [{
            "comparison_type": "capacity_n1",
            "before_analysis_id": steps[1]["analysis_id"],
            "after_analysis_id": steps[3]["analysis_id"],
            "before_baseline_max_congestion": 230.4,
            "after_baseline_max_congestion": 205.9,
            "relative_baseline_congestion_change": -0.1064,
            "counts": {
                "improved": 2,
                "regressed": 0,
                "unchanged": 6,
                "became_connected": 0,
                "became_disconnected": 0,
                "disconnected_both": 1,
            },
        }],
        "comparison_status": "completed",
        "workspace_modified": False,
        "answer": "Verified composite investigation completed.",
    }

    with sync_playwright() as playwright:
        browser = playwright.chromium.launch(headless=True)

        try:
            page = browser.new_page(
                viewport={"width": 1440, "height": 900}
            )
            page.goto(f"{BASE_URL}/ui")

            page.locator("#network").select_option(NETWORK_ID)

            page.route(
                f"**/api/v1/networks/{NETWORK_ID}/"
                "assistant/analyses/ask",
                lambda route: route.fulfill(
                    status=200,
                    content_type="application/json",
                    body=json.dumps(response_body),
                ),
            )

            page.locator("#jeylan-mode-experiment").click()
            page.locator(
                "#jeylan-experiment-question"
            ).fill(
                "Increase traffic, run N-1, upgrade edge 4, "
                "then run N-1 again."
            )

            page.locator("#jeylan-run-experiment").click(
                timeout=5000
            )

            panel = page.locator("#jeylan-composite-results")
            panel.wait_for(state="visible")

            assert page.locator(
                ".composite-step"
            ).count() == 4

            assert "−10.64%" not in panel.inner_text() or (
                    "10.64%" in panel.inner_text()
            )

            assert "2 congestion improvements" in (
                panel.inner_text()
            )

            page.locator(".composite-step").nth(1).click()

            assert steps[1]["analysis_id"] in (
                page.locator(
                    "#composite-step-detail"
                ).inner_text()
            )

            page.locator("#reset-network").click()

            assert panel.is_hidden()

        finally:
            browser.close()
