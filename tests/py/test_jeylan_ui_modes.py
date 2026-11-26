
"""Regression tests for Jeylan interaction modes and provenance."""

from __future__ import annotations

import os

from playwright.sync_api import sync_playwright


BASE_URL = os.environ.get(
    "ER_TEST_BASE_URL",
    "http://127.0.0.1:8000",
)


def test_jeylan_mode_switching():
    with sync_playwright() as playwright:
        browser = playwright.chromium.launch(headless=True)

        try:
            page = browser.new_page()
            page.goto(
                f"{BASE_URL}/ui",
                wait_until="domcontentloaded",
            )

            # Verify the default closed state.
            assert not page.locator("#jeylan-panel").is_visible()

            # Open Jeylan through the actual user interface.
            page.locator("#jeylan-launcher").click()

            ask_panel = page.locator("#jeylan-read-panel")
            experiment_panel = page.locator(
                "#jeylan-experiment-panel"
            )

            assert ask_panel.is_visible()
            assert not experiment_panel.is_visible()

            page.locator(
                "#jeylan-mode-experiment"
            ).click()

            assert experiment_panel.is_visible()
            assert not ask_panel.is_visible()

            question = page.locator(
                "#jeylan-experiment-question"
            )

            question.fill(
                "Increase capacity of edge 78 by 50%"
            )

            page.locator("#jeylan-mode-ask").click()

            assert ask_panel.is_visible()
            assert not experiment_panel.is_visible()

            # Switching tabs must preserve the draft.
            page.locator(
                "#jeylan-mode-experiment"
            ).click()

            assert question.input_value() == (
                "Increase capacity of edge 78 by 50%"
            )

        finally:
            browser.close()


def test_jeylan_provenance_and_workspace_reset():
    with sync_playwright() as playwright:
        browser = playwright.chromium.launch(headless=True)

        try:
            page = browser.new_page()
            page.goto(
                f"{BASE_URL}/ui",
                wait_until="domcontentloaded",
            )

            # Wait for the registered-network selector.
            page.wait_for_function("""
                () => {
                    const el = document.getElementById('network');
                    return el && !el.disabled && !!el.value;
                }
            """)

            # Open the detailed analysis results.
            page.locator(
                ".stitch-right-details > summary"
            ).click()

            provenance = page.locator(
                "#displayed-result-provenance"
            )

            assert not provenance.is_visible()

            page.evaluate("""
                setDisplayedResultProvenance({
                    origin: 'jeylan',
                    scenarioId: 'jeylan-test-001',
                    edgeId: 78,
                    capacityFactor: 1.5
                });
            """)

            assert provenance.is_visible()
            assert "Edge 78 ×1.5" in provenance.inner_text()
            assert "Workspace unchanged" in provenance.inner_text()

            # Open the parameter configuration.
            page.locator(
                "#stitch-parameter-disclosure > summary"
            ).click()

            # Open the workspace summary.
            page.locator(
                ".stitch-workspace-disclosure > summary"
            ).click()

            # Open the advanced network configuration.
            page.locator(
                ".stitch-native-advanced > summary"
            ).click()

            # Reset the workspace through the actual UI.
            page.locator("#reset-network").click()

            assert not provenance.is_visible()

        finally:
            browser.close()
