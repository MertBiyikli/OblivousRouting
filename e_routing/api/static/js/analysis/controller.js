function updateAnalysisMode() {
    clearDisplayedResultProvenance();
    networkGeneration++;


    const type = $('analysis-type').value;

    const capacity = type === 'capacity_upgrade';
    const failure = type === 'link_failure';
    const n1 = type === 'n1_resilience';
    const traffic = type === 'traffic_growth';

    const baselineDescriptions = {
        capacity_upgrade: 'Reference capacities · Current demand configuration',
        link_failure: 'Effective workspace · Before link failure',
        n1_resilience: 'Effective workspace · Before independent failures',
        traffic_growth: 'Effective capacities · Reference demand'
    };

    $('baseline-congestion-description').textContent = baselineDescriptions[type];

    // Show only configuration relevant to the selected mode.
    // Workspace configuration is independent of analysis selection.
    $('capacity-config').hidden = !capacity;
    $('failure-config').hidden = !failure;
    $('traffic-config').hidden = !traffic;
    $('execution-mode-config').hidden = !n1;

    // Demand model selection is available for all analysis types.
    $('demand').disabled = false;
    $('demand').closest('.field').hidden = false;

    // The shared demand multiplier is edited in traffic-config.
    // Keep the legacy scale input hidden to avoid two
    // conflicting controls for the same setting.
    $('scale').closest('.field').hidden = true;

    // Failure-specific visualization.
    $('failure-legend').hidden = !failure;


    // Invalidate results from the previous analysis mode.
    activeScenarioId = null;
    selectedLinkId = null;

    $('scenario-analysis').hidden = true;
    $('scenario-analysis').replaceChildren();

    $('reference-congestion').textContent = '—';
    $('upgraded-congestion').textContent = '—';
    $('capacity-improvement').textContent = '—';
    $('capacity-improvement').style.color = '';

    $('evidence').textContent = 'No analysis executed for this mode.';
    $('interpretation').textContent = 'Configure and execute an analysis to view results.';

    $('selected').textContent = 'Select a link after running an analysis.';

    $('ask-selected-link').disabled = true;
    $('jeylan-question').disabled = true;
    $('jeylan-ask').disabled = true;
    $('jeylan-status').textContent = 'Waiting for scenario execution…';

    $('status').textContent = 'Ready';
    clearError();

    setAnalysisMetricLabels(type);
    clearAnalysisSpecificMetrics();

    $('analysis-type-help').textContent = traffic ?
        'Reevaluate fixed routing under modified traffic demands.'
        : n1 ?
            'Evaluate independent physical-link failures.'
            : failure ?
                'Recompute routing after a physical-link failure.' :
                'Evaluate demand and capacity modifications.';

    $('run-help').textContent = traffic ?
        'Existing routing scheme · No optimization or topology modification.' :
        n1 ?
            'Selected demand model, seed 42. Recomputes routing for each link failure.' :
            failure ?
                'Recomputes modeled routing after a link failure.' :
                'Uses the registered persisted routing scheme.';

    $('run').textContent = traffic ?
        'Run traffic analysis' :
        n1 ?
            'Run N−1 analysis' :
            failure ?
                'Run failure analysis' :
                'Run simulation';

    invalidateJeylanContext();
}