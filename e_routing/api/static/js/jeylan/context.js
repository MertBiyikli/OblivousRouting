function activateJeylanContext(networkId, analysisId, analysisType) {
    activeAnalysisContext = {
        networkId,
        analysisId,
        analysisType
    };

    $('jeylan-question').disabled = false;
    $('jeylan-ask').disabled = false;
    $('ask-selected-link').disabled = lastLinks.length === 0;

    $('jeylan-status').textContent =
        'Analysis ready. Ask Jeylan about the results.';
}

function invalidateJeylanContext() {
    activeAnalysisContext = null;

    $('jeylan-question').disabled = true;
    $('jeylan-ask').disabled = true;
    $('ask-selected-link').disabled = true;

    $('jeylan-status').textContent =
        'Run an analysis before asking Jeylan.';
}