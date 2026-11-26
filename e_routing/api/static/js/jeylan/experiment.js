async function fetchJeylanAnalysisReport(
    networkId,
    analysisType,
    analysisId
) {
    const base =
        `/api/v1/networks/${encodeURIComponent(networkId)}`;

    let url;

    if (
        analysisType === 'capacity_upgrade' ||
        analysisType === 'link_failure'
    ) {
        url =
            `${base}/scenarios/` +
            `${encodeURIComponent(analysisId)}/analysis`;
    } else if (analysisType === 'n1_resilience') {
        url =
            `${base}/analyses/n1/` +
            encodeURIComponent(analysisId);
    } else if (analysisType === 'traffic_growth') {
        url =
            `${base}/analyses/traffic/` +
            encodeURIComponent(analysisId);
    } else {
        throw new Error('Unsupported Jeylan analysis type');
    }

    const response = await fetch(url);
    const body = await response.json();

    if (!response.ok) {
        throw new Error(
            typeof body.detail === 'string'
                ? body.detail
                : `Report unavailable (HTTP ${response.status})`
        );
    }

    if (body.status === 'running' || body.status === 'failed') {
        throw new Error(
            'The persisted analysis is not completed.'
        );
    }

    // Detailed reports may be enclosed in the API response.
    const report = body.analysis ?? body;

    const valid =
        analysisType === 'n1_resilience'
            ? report.analysis_kind === 'n1_resilience'
            : report.scenario_kind === analysisType;

    if (!valid) {
        throw new Error('Persisted analysis type mismatch');
    }

    return report;
}


async function displayJeylanAnalysisResult({
                                               networkId,
                                               analysisId,
                                               analysisType,
                                               report,
                                               stillCurrent
                                           }) {
    if (!stillCurrent()) return;

    selectedLinkId = null;
    activeScenarioId = null;
    setAnalysisMetricLabels(analysisType);
    clearAnalysisSpecificMetrics();
    $('baseline-congestion-description').textContent = ({
        capacity_upgrade: 'Reference capacities · Current demand configuration',
        link_failure: 'Effective workspace · Before link failure',
        n1_resilience: 'Effective workspace · Before independent failures',
        traffic_growth: 'Effective capacities · Reference demand'
    })[analysisType];
    $('congestion-delta').style.color = '';

    if (analysisType === 'capacity_upgrade') {
        await displayCapacityAnalysisResult({
            networkId,
            scenarioId: analysisId,
            report,
            updateInterpretation: false,
            stillCurrent
        });

    } else if (analysisType === 'traffic_growth') {
        if (
            !Array.isArray(report.baseline?.links) ||
            !Array.isArray(report.modified?.links)
        ) {
            throw new Error('Traffic report has no link measurements');
        }

        renderTrafficReport(report);

    } else if (analysisType === 'n1_resilience') {
        if (
            !Array.isArray(report.baseline?.links) ||
            !Array.isArray(report.outcomes)
        ) {
            throw new Error('N−1 report is incomplete');
        }

        const baseline = report.baseline;
        lastLinks = baseline.links;

        $('baseline-congestion').textContent =
            Number(baseline.max_congestion).toFixed(3);

        $('scenario-congestion').textContent = 'N/A';
        $('congestion-delta').textContent = 'N/A';
        $('congestion').textContent =
            Number(baseline.max_congestion).toFixed(3);

        $('count').textContent = String(lastLinks.length);
        $('overloaded').textContent = String(
            lastLinks.filter(link => link.utilization > 1).length
        );

        $('failure-legend').hidden = true;
        $('topology').textContent =
            'Network topology · N−1 baseline';

        drawGraph(lastLinks);
        renderBottleneckRanking(lastLinks);

        renderN1Analysis(report, outcome => {
            if (!stillCurrent()) return;

            selectedLinkId = outcome.edge_id;

            const successful =
                outcome.status === 'success' &&
                Array.isArray(outcome.reoptimized?.links);

            const links = successful
                ? outcome.reoptimized.links
                : baseline.links.filter(
                    link => link.edge_id !== outcome.edge_id
                );

            lastLinks = links;

            drawGraph(links, {
                baselineLinks: baseline.links,
                failedIds: [outcome.edge_id]
            });

            renderBottleneckRanking(links);
            $('failure-legend').hidden = false;

            $('topology').textContent =
                `Network topology · Edge ${outcome.edge_id} failed`;

            // For disconnects and execution failures these are
            // baseline link colors only, not post-failure measurements.
            $('count').textContent = successful ? String(links.length) : 'N/A';
            $('overloaded').textContent = successful
                ? String(links.filter(link => Number(link.utilization) > 1).length)
                : 'N/A';
            $('reference-congestion').textContent = `Edge ${outcome.edge_id}`;
            $('upgraded-congestion').textContent = successful
                ? Number(outcome.reoptimized.max_congestion).toFixed(3)
                : 'N/A';
            $('capacity-improvement').textContent =
                outcome.status === 'disconnected' ? 'Disconnected' :
                    successful && Number(baseline.max_congestion) > 0
                        ? `${((outcome.reoptimized.max_congestion /
                            baseline.max_congestion - 1) * 100).toFixed(2)}%`
                        : 'Execution failed';
            $('third-help-1').textContent =
                `${outcome.source} → ${outcome.target} · removed link`;
            $('third-help-2').textContent = successful
                ? 'Recomputed maximum congestion'
                : 'No post-failure congestion is defined';
            $('third-help-3').textContent =
                'Relative congestion change or failure status';

            if (successful) {
                const after = outcome.reoptimized.max_congestion;
                const before = baseline.max_congestion;

                $('scenario-congestion').textContent =
                    Number(after).toFixed(3);
                $('congestion').textContent =
                    Number(after).toFixed(3);

                $('congestion-delta').textContent =
                    before > 0
                        ? `${((after / before - 1) * 100).toFixed(2)}%`
                        : 'N/A';

                $('selected').textContent =
                    `Edge ${outcome.edge_id}: routing recomputed.`;
            } else {
                $('scenario-congestion').textContent = 'N/A';
                $('congestion').textContent = 'N/A';
                $('congestion-delta').textContent = 'N/A';

                $('selected').textContent =
                    outcome.status === 'disconnected'
                        ? `Edge ${outcome.edge_id}: network disconnected. ` +
                        'Surviving links show baseline reference values.'
                        : `Edge ${outcome.edge_id}: execution failed. ` +
                        'No post-failure measurements available.';
            }
        });

    } else if (analysisType === 'link_failure') {
        const solver = report.solver_results?.[0];

        if (!solver?.baseline) {
            throw new Error('Missing link-failure baseline');
        }

        const baseline = solver.baseline;
        const modified = solver.reoptimized;
        const failedId = report.link_failures?.[0];

        if (!Number.isInteger(failedId)) {
            throw new Error('Missing failed-link identity');
        }

        const successful = Array.isArray(modified?.links);

        lastLinks = successful
            ? modified.links
            : baseline.links.filter(
                link => link.edge_id !== failedId
            );

        selectedLinkId = failedId;

        drawGraph(lastLinks, {
            baselineLinks: baseline.links,
            failedIds: [failedId]
        });

        renderBottleneckRanking(lastLinks);
        $('failure-legend').hidden = false;

        $('topology').textContent =
            `Network topology · Edge ${failedId} failed`;

        $('baseline-congestion').textContent =
            Number(baseline.max_congestion).toFixed(3);

        $('scenario-congestion').textContent = successful
            ? Number(modified.max_congestion).toFixed(3)
            : 'N/A';

        $('congestion').textContent = successful
            ? Number(modified.max_congestion).toFixed(3)
            : 'N/A';

        $('congestion-delta').textContent =
            successful && baseline.max_congestion > 0
                ? (
                (modified.max_congestion /
                    baseline.max_congestion - 1) * 100
            ).toFixed(2) + '%'
                : 'N/A';

        $('count').textContent = String(lastLinks.length);
        $('overloaded').textContent = successful
            ? String(
                lastLinks.filter(link => link.utilization > 1).length
            )
            : 'N/A';

        $('selected').textContent = successful
            ? `Link ${failedId} removed. Routing recomputed.`
            : `Link ${failedId} disconnects the network. ` +
            'No post-failure congestion is defined.';

        const failedLink = baseline.links.find(link => link.edge_id === failedId);
        const criticalLink = successful
            ? [...lastLinks].sort((a, b) => b.utilization - a.utilization)[0]
            : null;
        const largestShift = successful
            ? [...(solver.link_deltas || [])].sort((a, b) =>
                Math.abs(b.delta_from_baseline) -
                Math.abs(a.delta_from_baseline))[0]
            : null;

        $('reference-congestion').textContent = `Edge ${failedId}`;
        $('upgraded-congestion').textContent = criticalLink
            ? `Edge ${criticalLink.edge_id}` : 'N/A';
        $('capacity-improvement').textContent = largestShift
            ? `${largestShift.delta_from_baseline >= 0 ? '+' : ''}` +
            `${Number(largestShift.delta_from_baseline).toFixed(3)}`
            : 'N/A';
        $('third-help-1').textContent = failedLink
            ? `${failedLink.source} → ${failedLink.target} · out of service`
            : 'Removed physical link';
        $('third-help-2').textContent = criticalLink
            ? `${criticalLink.source} → ${criticalLink.target} · ` +
            `${Number(criticalLink.utilization).toFixed(3)}× capacity`
            : 'No recomputed routing metrics';
        $('third-help-3').textContent = largestShift
            ? `Edge ${largestShift.edge_id} · absolute load change`
            : 'No post-failure load deltas';

        await renderScenarioAnalysis(networkId, analysisId);
    }

    if (!stillCurrent()) return;

    $('status').textContent = 'Completed';
    $('evidence').textContent =
        `Jeylan simulated analysis · ${analysisType} · ` +
        `Network: ${networkId} · Report: ${analysisId}`;

    activateJeylanContext(
        networkId,
        analysisId,
        analysisType
    );
}