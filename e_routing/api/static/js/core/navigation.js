function switchStudioView(view) {
    topologyActive = view === 'topology';
    if (topologyActive && window.stitchRefreshTopology) requestAnimationFrame(window.stitchRefreshTopology);
    document.body.classList.toggle('topology-view', topologyActive);
    // Keep the Analysis graph outside the topology comparison DOM entirely.
    // CSS-only hiding left a two-pane SVG layout affecting the analysis viewport.
    const graphSvg = $('graph');
    const analysisGraphWrap = document.querySelector('.graph-panel .graph-wrap');
    const comparisonCurrentPane = $('stitch-current-pane');
    if (graphSvg && analysisGraphWrap && comparisonCurrentPane) {
        if (topologyActive) {
            if (graphSvg.parentElement !== comparisonCurrentPane) comparisonCurrentPane.appendChild(graphSvg);
        } else {
            if (graphSvg.parentElement !== analysisGraphWrap) analysisGraphWrap.appendChild(graphSvg);
            const baselineSvg = $('stitch-baseline-graph');
            if (baselineSvg) baselineSvg.replaceChildren();
        }
    }

    for (const name of ['topology', 'analysis']) {
        const active = name === view;
        $('nav-' + name).classList.toggle('active', active);
        $('nav-' + name).setAttribute('aria-pressed', String(active));
    }
    if (topologyActive) updateTopologyInspector();
}