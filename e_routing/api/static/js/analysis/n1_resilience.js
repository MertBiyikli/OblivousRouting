async function executeN1UI() {
    clearDisplayedResultProvenance();
    clearError();

    const networkId = $('network').value;
    const executionMode = $('execution-mode').value;
    if (!networkId || $('network').disabled) {
        showError('Select a registered network.');
        return;
    }

    const generation = networkGeneration;
    const stillCurrent = () => generation === networkGeneration && $('network').value === networkId && $('analysis-type').value === 'n1_resilience';

    const analysisId = 'n1-' + (globalThis.crypto?.randomUUID?.() || String(Date.now()) + Math.random().toString(36).slice(2))
        .replace(/[^A-Za-z0-9_-]/g, '');

    const button = $('run');
    button.disabled = true;
    button.textContent = 'Evaluating all links…';

    $('status').textContent = 'Running';
    $('evidence').textContent = 'Running N−1 resilience analysis…';
    $('scenario-analysis').hidden = true;
    $('scenario-analysis').replaceChildren();

    $('graph').replaceChildren();
    $('empty').hidden = false;
    $('empty').textContent = 'Analyzing independent link failures…';
    $('large-network-analysis').hidden = true;
    $('failure-legend').hidden = true;

    // N-1 reports are not normal ScenarioRepository/Jeylan scenarios.
    activeScenarioId = null;
    selectedLinkId = null;
    lastLinks = [];
    $('jeylan-question').disabled = true;
    $('jeylan-ask').disabled = true;
    $('ask-selected-link').disabled = true;
    $('jeylan-status').textContent = 'N−1 interpretation is not connected to Jeylan yet.';

    for (const id of ['congestion', 'overloaded', 'count', 'baseline-congestion', 'scenario-congestion', 'congestion-delta', 'reference-congestion', 'upgraded-congestion', 'capacity-improvement']) {
        $(id).textContent = '—';
    }

    try {
        const base = `/api/v1/networks/${encodeURIComponent(networkId)}` + `/analyses/n1`;

        const response = await fetch(base, {
            method: 'POST', headers: {'Content-Type': 'application/json'}, body: JSON.stringify({
                analysis_id: analysisId,
                execution_mode: executionMode,
                demand_model: $('demand').value,
                seed: 42,
                workspace: currentWorkspaceSnapshot()
            })
        });

        const result = await response.json();

        if (!response.ok) {
            throw new Error(result.detail || `N−1 execution failed (HTTP ${response.status})`);
        }

        if (result.status !== 'completed' || result.analysis_id !== analysisId) {
            throw new Error('Incomplete N−1 execution response');
        }

        if (!stillCurrent()) return;

        const detailResponse = await fetch(`${base}/${encodeURIComponent(analysisId)}`);

        if (!detailResponse.ok) {
            throw new Error(`N−1 report unavailable (HTTP ${detailResponse.status})`);
        }

        const report = await detailResponse.json();

        if (!stillCurrent()) return;

        if (report.analysis_kind !== 'n1_resilience' || !Array.isArray(report.outcomes) || !Array.isArray(report.baseline?.links)) {
            throw new Error('Invalid N−1 analysis report');
        }


        setAnalysisMetricLabels('n1_resilience');
        const baseline = report.baseline;
        const before = baseline.max_congestion;

        $('baseline-congestion').textContent = fmt(before, 3);
        $('status').textContent = 'Completed';

        $('evidence').textContent = `Network: ${networkId} · N−1: ${analysisId} · ` + `Mode: ${executionMode.toUpperCase()} · ` + `Solver: ${report.solver} · ` + `Fingerprint: ${report.graph_fingerprint}`;

        $('interpretation').textContent = `${report.summary.tested_links} independent physical-link ` + `failures evaluated. Select a failure to inspect its ` + `network-level consequences.`;


        function selectFailure(outcome) {
            // Do not reuse a previous scenario's Jeylan context.
            selectedLinkId = outcome.edge_id;
            $('ask-selected-link').disabled = !activeAnalysisContext;

            const disconnected = outcome.status === 'disconnected';
            const successful = outcome.status === 'success';
            const modified = successful ? outcome.reoptimized : null;

            // If disconnected, show baseline links only as topology
            // context—not as a fictional post-failure routing result.
            const links = successful ? modified.links.map(link => {
                const previous = baseline.links.find(item => item.edge_id === link.edge_id);
                return {
                    ...link,
                    baseline_utilization: previous?.utilization ?? null,
                    baseline_capacity: previous?.capacity ?? null,
                    utilization_change: previous ? link.utilization - previous.utilization : null
                };
            }) : baseline.links.filter(link => link.edge_id !== outcome.edge_id);

            lastLinks = links;

            drawGraph(links, {
                baselineLinks: baseline.links, failedIds: [outcome.edge_id]
            });

            $('failure-legend').hidden = false;
            renderBottleneckRanking(links);

            $('topology').textContent = `Network topology · Edge ${outcome.edge_id} failed`;

            $('selected').textContent = disconnected ? `Edge ${outcome.edge_id} disconnects demanded endpoints. ` + `Surviving links show baseline reference values only.` : successful ? `Edge ${outcome.edge_id} removed. ` + `Surviving links show recomputed utilization.` : `Edge ${outcome.edge_id}: recomputation failed. ` + `No post-failure metrics available.`;

            const after = modified?.max_congestion ?? null;
            const change = successful && before > 0 ? 100 * (after - before) / before : null;

            $('scenario-congestion').textContent = fmt(after, 3);
            $('congestion').textContent = fmt(after, 3);
            $('congestion-delta').textContent = change === null ? 'N/A' : `${change >= 0 ? '+' : ''}${fmt(change)}%`;

            $('congestion-delta').style.color = change === null ? '#94a3b8' : change < 0 ? '#22c55e' : change > 0 ? '#ef4444' : '#94a3b8';

            $('count').textContent = successful ? String(links.length) : 'N/A';
            $('overloaded').textContent = successful ? String(links.filter(l => l.utilization > 1).length) : 'N/A';

            $('reference-congestion').textContent = `Edge ${outcome.edge_id}`;
            $('upgraded-congestion').textContent = successful ? fmt(after, 3) : 'N/A';
            $('capacity-improvement').textContent = disconnected ? 'Disconnected' : successful ? (change === null ? 'N/A' : `${change >= 0 ? '+' : ''}${fmt(change)}%`) : 'Execution failed';

            $('third-help-1').textContent = `${outcome.source} → ${outcome.target} · removed link`;
            $('third-help-2').textContent = 'Maximum post-failure congestion';
            $('third-help-3').textContent = 'Relative congestion change or failure status';
        }

        renderN1Analysis(report, selectFailure);
        activateJeylanContext(networkId, analysisId, 'n1_resilience');

    } catch (error) {
        if (stillCurrent()) {
            $('status').textContent = 'Error';
            $('evidence').textContent = 'N−1 analysis did not complete.';
            showError(error.message || String(error));
        }
    } finally {
        if (stillCurrent()) {
            button.disabled = false;
            button.textContent = 'Run N−1 analysis';
        }
    }
}