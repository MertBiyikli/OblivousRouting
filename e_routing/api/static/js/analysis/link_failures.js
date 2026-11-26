async function executeLinkFailureUI() {
    clearDisplayedResultProvenance();
    clearError();

    const networkId = $('network').value;
    const rawEdge = $('failed-edge').value.trim();
    const edgeId = Number(rawEdge);

    if (!networkId || $('network').disabled) {
        showError('Select a registered network.');
        return;
    }

    if (rawEdge === '' || !Number.isSafeInteger(edgeId) || edgeId < 0) {
        showError('Select a valid physical link ID.');
        return;
    }

    const scenarioId = 'failure-' + (globalThis.crypto?.randomUUID?.() || String(Date.now()) + Math.random().toString(36).slice(2)).replace(/[^A-Za-z0-9_-]/g, '');
    const generation = networkGeneration;

    const stillCurrent = () => generation === networkGeneration && $('network').value === networkId && $('analysis-type').value === 'link_failure';

    const request = {
        scenario_id: scenarioId,
        scenario_kind: 'link_failure',
        demand_model: $('demand').value,
        workspace: currentWorkspaceSnapshot(),
        link_failures: [edgeId]
    };

    const button = $('run');
    button.disabled = true;
    button.textContent = 'Analyzing failure…';

    $('status').textContent = 'Running';
    $('scenario-analysis').hidden = true;
    $('scenario-analysis').replaceChildren();

    // The prior scenario must not remain active.
    activeScenarioId = null;
    selectedLinkId = null;
    lastLinks = [];

    $('jeylan-question').disabled = true;
    $('jeylan-ask').disabled = true;
    $('ask-selected-link').disabled = true;
    $('jeylan-status').textContent = 'Detailed failure analysis is not connected to Jeylan yet.';

    $('graph').replaceChildren();
    $('empty').hidden = false;
    $('empty').textContent = 'Running link failure analysis…';

    $('large-network-analysis').hidden = true;
    $('bottleneck-list').replaceChildren();

    for (const id of ['congestion', 'overloaded', 'count', 'baseline-congestion', 'scenario-congestion', 'congestion-delta', 'reference-congestion', 'upgraded-congestion', 'capacity-improvement']) {
        $(id).textContent = '—';
    }

    $('evidence').textContent = 'Executing failure scenario…';
    $('interpretation').textContent = 'Recomputing routing on the modified network.';

    try {
        const url = `/api/v1/networks/${encodeURIComponent(networkId)}` + `/scenarios/detailed-analysis`;

        const response = await fetch(url, {
            method: 'POST', headers: {'Content-Type': 'application/json'}, body: JSON.stringify(request)
        });

        const body = await response.json();

        if (!stillCurrent()) return;

        if (response.ok && body.status === 'disconnected') {
            const failedId = Number(body.failed_edge_id);
            const baselineLinks = baselineTopologyCache.get(networkId);

            $('status').textContent = 'Disconnected';
            $('topology').textContent = `Network topology · Link ${failedId} disconnected`;

            $('failure-legend').hidden = false;

            if (baselineLinks?.length) {
                // No valid post-failure utilization exists.
                // Draw surviving links as unmeasured.
                const survivingLinks = baselineLinks
                    .filter(link => link.edge_id !== failedId)
                    .map(link => ({
                        ...link, utilization: NaN, load: null
                    }));

                drawGraph(survivingLinks, {
                    baselineLinks, failedIds: [failedId]
                });

                $('count').textContent = String(survivingLinks.length);

                $('selected').textContent = `Dashed amber: failed edge ${failedId}. ` + 'Surviving links have no post-failure routing measurements.';
            } else {
                $('graph').replaceChildren();
                $('empty').hidden = false;
                $('empty').textContent = `Network disconnected after failure of link ${failedId}. ` + 'Original topology has not yet been loaded.';
                $('failure-legend').hidden = true;
                $('selected').textContent = 'Run a successful scenario on this network to load its topology.';
            }

            $('interpretation').textContent = `Physical link ${failedId} is critical: ` + 'its removal disconnects endpoints with positive traffic demand. ' + 'Post-failure congestion is not defined.';

            $('evidence').textContent = `Scenario: ${body.scenario_id} · ` + `Network: ${networkId} · Disconnected topology`;

            $('scenario-analysis').hidden = false;
            $('scenario-analysis').textContent = body.reason || 'The failed link disconnects demanded endpoints.';

            return;
        }

        if (!response.ok) {
            const detail = body.detail;
            throw new Error(typeof detail === 'string' ? detail : JSON.stringify(detail || body));
        }

        if (!stillCurrent()) return;

        if (body.status !== 'completed' || !body.report_available || body.analysis?.scenario_kind !== 'link_failure') {
            throw new Error('Incomplete link failure response');
        }

        setAnalysisMetricLabels('link_failure');
        const report = body.analysis;
        const solver = report.solver_results?.[0];

        if (!solver) {
            throw new Error('Missing solver analysis results');
        }

        const baseline = solver.baseline;
        const modified = solver.reoptimized;

        if (Array.isArray(baseline.links) && baseline.links.length > 0) {
            baselineTopologyCache.set(networkId, baseline.links.map(link => ({...link})));
        }

        const before = baseline.max_congestion;
        const after = modified.max_congestion;
        const change = before > 0 ? 100 * (after - before) / before : null;

        const format = value => Number.isFinite(value) ? value.toFixed(3) : '—';

        $('baseline-congestion').textContent = format(before);
        $('scenario-congestion').textContent = format(after);
        $('congestion').textContent = format(after);

        $('congestion-delta').textContent = change === null ? 'N/A' : `${change > 0 ? '+' : ''}${change.toFixed(2)}%`;

        $('congestion-delta').style.color = change === null ? '#94a3b8' : change < 0 ? '#22c55e' : change > 0 ? '#ef4444' : '#94a3b8';

        const baselineByEdge = new Map(baseline.links.map(link => [link.edge_id, link]));
        lastLinks = modified.links.map(link => {
            const previous = baselineByEdge.get(link.edge_id);
            return {
                ...link,
                baseline_utilization: previous?.utilization ?? null,
                baseline_capacity: previous?.capacity ?? null,
                utilization_change: previous ? link.utilization - previous.utilization : null
            };
        });
        const links = lastLinks;
        $('failure-legend').hidden = false;
        drawGraph(links, {
            baselineLinks: baseline.links, failedIds: report.link_failures
        });
        renderBottleneckRanking(links);
        $('topology').textContent = `Network topology · Link ${edgeId} failed`;
        $('selected').textContent = `Dashed amber: failed edge ${edgeId}. Click an edge to inspect it.`;

        $('count').textContent = String(links.length);
        $('overloaded').textContent = String(links.filter(link => link.utilization > 1).length);

        // The third metric row is analysis-type-specific.
        const failedLink = baselineByEdge.get(edgeId);
        const bottleneck = [...links].sort((a, b) => b.utilization - a.utilization)[0];
        const largestShift = [...(solver.link_deltas || [])].sort((a, b) => Math.abs(b.delta_from_baseline) - Math.abs(a.delta_from_baseline))[0];
        $('reference-congestion').textContent = failedLink ? `Edge ${edgeId}` : 'N/A';
        $('upgraded-congestion').textContent = bottleneck ? `Edge ${bottleneck.edge_id}` : 'N/A';
        $('capacity-improvement').textContent = largestShift ? `${largestShift.delta_from_baseline >= 0 ? '+' : ''}` + `${format(largestShift.delta_from_baseline)}` : 'N/A';
        $('capacity-improvement').style.color = '';
        $('third-help-1').textContent = failedLink ? `${failedLink.source} → ${failedLink.target} · out of service` : 'Failed link unavailable in baseline';
        $('third-help-2').textContent = bottleneck ? `${bottleneck.source} → ${bottleneck.target} · ${format(bottleneck.utilization)}× capacity` : 'No surviving links';
        $('third-help-3').textContent = largestShift ? `Edge ${largestShift.edge_id} · absolute load change` : 'No surviving link deltas';

        $('status').textContent = 'Completed';

        $('evidence').textContent = `Scenario: ${scenarioId} · ` + `Network: ${networkId} · ` + `Solver: ${solver.solver_algorithms.join(', ')} · ` + `Mode: ${report.execution_mode}`;

        $('interpretation').textContent = `Physical link ${edgeId} was removed. ` + `Maximum congestion changed from ${format(before)} ` + `to ${format(after)}. ` + `This is a modeled rerouting result.`;

        $('scenario-analysis').hidden = false;

        await renderScenarioAnalysis(networkId, body.scenario_id);

        if (!stillCurrent()) return;

        activateJeylanContext(networkId, scenarioId, 'link_failure');


        if (!stillCurrent()) {
            $('scenario-analysis').hidden = true;
            $('scenario-analysis').replaceChildren();
        }
    } catch (error) {
        if (stillCurrent()) {
            $('status').textContent = 'Error';
            $('evidence').textContent = 'Failure analysis did not complete.';
            showError(error.message || String(error));
        }
    } finally {
        if (stillCurrent()) {
            button.disabled = false;
            button.textContent = 'Run failure analysis';
        }
    }
}