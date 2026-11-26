async function displayCapacityAnalysisResult({
                                                 networkId,
                                                 scenarioId,
                                                 report,
                                                 updateInterpretation = true,
                                                 stillCurrent = () => true
                                             }) {
    if (report?.scenario_kind !== 'capacity_upgrade'
        || !Array.isArray(report.solver_results)
        || report.solver_results.length !== 1) {
        throw new Error('Invalid capacity analysis report');
    }

    setAnalysisMetricLabels('capacity_upgrade');
    const solver = report.solver_results[0];

    if (!solver?.baseline || !solver?.fixed_routing || !solver?.reoptimized || !Array.isArray(solver.baseline.links) || !Array.isArray(solver.reoptimized.links)) {
        throw new Error('Missing three-way capacity comparison');
    }
    /*
    for(let i = 0; i < report.solver_results.length; i++) {
        const solver = report.solver_results[i];
        if (!solver?.baseline || !solver?.fixed_routing || !solver?.reoptimized || !Array.isArray(solver.baseline.links) || !Array.isArray(solver.reoptimized.links)) {
            throw new Error('Missing three-way capacity comparison');
        }
    }*/

    if (!stillCurrent()) return;

    const baseline = solver.baseline;
    const fixed = solver.fixed_routing;
    const optimized = solver.reoptimized;

    const before = baseline.max_congestion;
    const fixedCongestion = fixed.max_congestion;
    const after = optimized.max_congestion;

    const format = value => Number.isFinite(value) ? value.toFixed(3) : '—';

    const change = before > 0 ? 100 * (after - before) / before : null;


    // First metric rows.
    $('baseline-congestion').textContent = format(before);
    $('scenario-congestion').textContent = format(after);
    $('congestion').textContent = format(after);

    $('baseline-congestion-description').textContent = 'Reference capacities · Current demand configuration';

    $('congestion-delta').textContent = change === null ? 'N/A' : `${change >= 0 ? '+' : ''}${change.toFixed(2)}%`;

    $('congestion-delta').style.color = change === null ? '#94a3b8' : change < 0 ? '#22c55e' : change > 0 ? '#ef4444' : '#94a3b8';

    // Three-way capacity comparison.
    $('reference-congestion').textContent = format(before);
    $('upgraded-congestion').textContent = format(fixedCongestion);
    $('capacity-improvement').textContent = format(after);
    $('capacity-improvement').style.color = '';

    $('third-label-1').textContent = 'Baseline routing';
    $('third-label-2').textContent = 'Capacity-only';
    $('third-label-3').textContent = 'Reoptimized routing';

    $('third-help-1').textContent = 'Original capacities';
    $('third-help-2').textContent = 'Fixed routing · Upgraded capacities';
    $('third-help-3').textContent = 'Routing recomputed after capacity upgrade';

    // Build link-level differences for topology visualization.
    const baselineByEdge = new Map(baseline.links.map(link => [link.edge_id, link]));

    lastLinks = optimized.links.map(link => {
        const previous = baselineByEdge.get(link.edge_id);

        return {
            ...link,
            baseline_utilization: previous?.utilization ?? null,
            utilization_change: previous ? link.utilization - previous.utilization : null
        };
    });

    selectedLinkId = null;
    $('selected').textContent = 'Select a link to inspect its utilization.';

    drawGraph(lastLinks);
    renderBottleneckRanking(lastLinks);

    $('topology').textContent = 'Network topology · Reoptimized Capacity Planning';

    $('failure-legend').hidden = true;
    $('count').textContent = String(lastLinks.length);

    $('overloaded').textContent = String(lastLinks.filter(link => link.utilization > 1).length);

    $('status').textContent = 'Completed';

    $('evidence').textContent = `Scenario: ${scenarioId} · ` + `Network: ${networkId} · ` + `Solver: ${solver.solver_algorithms.join(', ')} · ` + `Mode: ${report.execution_mode}`;

    if (updateInterpretation) {
        $('interpretation').textContent = `Capacity planning completed. Maximum congestion ` + `changed from ${format(before)} to ${format(after)} ` + `after routing reoptimization ` + `(${change === null ? 'N/A' : change.toFixed(2) + '%'}). ` + `The capacity-only result was ${format(fixedCongestion)}.`;
    }

    // Reuse the existing detailed-report visualization.
    const panel = $('scenario-analysis');
    panel.hidden = false;

    await renderScenarioAnalysis(networkId, scenarioId);

    if (!stillCurrent()) {
        panel.hidden = true;
        panel.replaceChildren();
        return;
    }

    // Enable read-only Jeylan questions for the displayed experiment.
    activateJeylanContext(networkId, scenarioId, 'capacity_upgrade');

    setDisplayedResultProvenance({
        origin: 'studio', scenarioId
    });
}