from pathlib import Path


UI_DIRECTORY = Path(__file__).resolve().parent / "static"
TEMPLATE_DIRECTORY = UI_DIRECTORY / "templates"

UI_FRAGMENTS = {
    "<!-- INCLUDE:HEADER -->": "header.html",
    "<!-- INCLUDE:ANALYSIS_SIDEBAR -->": "analysis-sidebar.html",
    "<!-- INCLUDE:ANALYSIS_WORKSPACE -->": "analysis-workspace.html",
    "<!-- INCLUDE:JEYLAN_PANEL -->": "jeylan-panel.html",
    "<!-- INCLUDE:HISTORY_PAGE -->": "history-page.html",
    "<!-- INCLUDE:JEYLAN_LAUNCHER -->": "jeylan-launcher.html",
}


def render_ui() -> str:
    html = (UI_DIRECTORY / "index.html").read_text(
        encoding="utf-8"
    )

    for marker, filename in UI_FRAGMENTS.items():
        if html.count(marker) != 1:
            raise RuntimeError(
                f"Expected exactly one UI marker: {marker}"
            )

        fragment = (
                TEMPLATE_DIRECTORY / filename
        ).read_text(encoding="utf-8")

        html = html.replace(marker, fragment, 1)

    return html