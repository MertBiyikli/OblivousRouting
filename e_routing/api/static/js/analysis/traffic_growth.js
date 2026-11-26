function renderTrafficOverrides() {
    const container = $('traffic-overrides');
    container.replaceChildren();

    if (!trafficOverrides.size) {
        container.textContent = 'Commodity overrides: none';
        return;
    }

    for (const [key, factor] of trafficOverrides) {
        const row = document.createElement('div');
        const label = document.createElement('span');
        const remove = document.createElement('button');

        label.textContent = `${key.replace(':', ' → ')} ×${factor}`;

        remove.type = 'button';
        remove.className = 'secondary';
        remove.textContent = 'Remove';

        remove.addEventListener('click', () => {
            trafficOverrides.delete(key);
            renderTrafficOverrides();
            refreshWorkspaceSummary();
        });

        row.append(label, remove);
        container.appendChild(row);
    }
}

function trafficPercent(value) {
    return Number.isFinite(value) ? `${value >= 0 ? '+' : ''}${(100 * value).toFixed(2)}%` : 'N/A';
}

function renderTrafficReport(report) {
    setAnalysisMetricLabels('traffic_growth');
    const baseline = report.baseline;
    const modified = report.modified;

    const baselineLinks = new Map(baseline.links.map(link => [link.edge_id, link]));

    lastLinks = modified.links.map(link => {
        const previous = baselineLinks.get(link.edge_id);

        return {
            ...link,
            baseline_load: previous?.load ?? null,
            baseline_utilization: previous?.utilization ?? null,
            baseline_capacity: previous?.capacity ?? null,
            utilization_change: previous ? link.utilization - previous.utilization : null
        };
    });

    selectedLinkId = null;

    $('selected').textContent = 'Traffic analysis updated. Select a link to inspect ' + 'its baseline and modified load.';
    $('topology').textContent = 'Network topology · Traffic Growth';
    $('failure-legend').hidden = true;
    $('baseline-congestion').textContent = baseline.max_congestion.toFixed(3);
    $('scenario-congestion').textContent = modified.max_congestion.toFixed(3);
    $('congestion-delta').textContent = trafficPercent(report.relative_congestion_change);
    $('reference-congestion').textContent = baseline.total_demand_volume.toLocaleString(undefined, {maximumFractionDigits: 2});
    $('upgraded-congestion').textContent = modified.total_demand_volume.toLocaleString(undefined, {maximumFractionDigits: 2});
    $('capacity-improvement').textContent = trafficPercent(report.relative_demand_change);
    $('congestion').textContent = modified.max_congestion.toFixed(3);
    $('overloaded').textContent = String(lastLinks.filter(link => link.utilization > 1).length);
    $('count').textContent = String(lastLinks.length);

    drawGraph(lastLinks);
    renderBottleneckRanking(lastLinks);

    const panel = $('scenario-analysis');
    panel.hidden = false;
    panel.replaceChildren();

    const heading = document.createElement('h2');
    heading.textContent = 'Traffic Growth Analysis';
    panel.appendChild(heading);

    const summary = [`Baseline demand volume: ${baseline.total_demand_volume.toFixed(2)}`, `Modified demand volume: ${modified.total_demand_volume.toFixed(2)}`, `Demand change: ${trafficPercent(report.relative_demand_change)}`, `Maximum congestion change: ${trafficPercent(report.relative_congestion_change)}`, `Global multiplier: ${report.requested_global_scale}`, `Commodity overrides: ${report.commodity_overrides.length}`];

    for (const line of summary) {
        const paragraph = document.createElement('p');
        paragraph.textContent = line;
        panel.appendChild(paragraph);
    }


    if (report.commodity_overrides.length > 0) {
        const heading = document.createElement('h3');
        heading.textContent = 'Commodity transformations';
        panel.appendChild(heading);

        const list = document.createElement('ul');

        for (const item of report.commodity_overrides) {
            const row = document.createElement('li');

            row.textContent = `${item.source} → ${item.target}: ` + `${Number(item.factor).toFixed(2)}× demand`;

            list.appendChild(row);
        }

        panel.appendChild(list);
    }


    const headingLinks = document.createElement('h3');
    headingLinks.textContent = 'Largest utilization changes';
    panel.appendChild(headingLinks);

    const ranked = [...lastLinks]
        .sort((a, b) => Math.abs(b.utilization_change ?? 0) - Math.abs(a.utilization_change ?? 0))
        .slice(0, 10);

    for (const link of ranked) {
        const paragraph = document.createElement('p');
        paragraph.textContent = `Edge ${link.edge_id} (${link.source} → ${link.target}): ` + `${link.baseline_utilization.toFixed(3)} → ` + `${link.utilization.toFixed(3)} ` + `(Δ ${(link.utilization_change >= 0 ? '+' : '')}` + `${link.utilization_change.toFixed(3)})`;

        panel.appendChild(paragraph);
    }

    $('status').textContent = 'Completed';
    $('interpretation').textContent = `Fixed-routing analysis completed. Total traffic changed by ` + `${trafficPercent(report.relative_demand_change)}, while maximum ` + `congestion changed by ${trafficPercent(report.relative_congestion_change)}. ` + `No routing optimization or network modification was performed.`;
}

async function executeTrafficGrowthUI() {
    clearDisplayedResultProvenance();
    clearError();

    const networkId = $('network').value;

    if (!networkId || $('network').disabled) {
        showError('Select a registered network.');
        return;
    }

    const scale = Number($('traffic-scale').value);

    if ($('traffic-scale').value.trim() === '' || !Number.isFinite(scale) || scale <= 0 || scale > 100) {
        showError('Enter a global multiplier greater than zero and at most 100.');
        return;
    }

    const executionGeneration = networkGeneration;

    const stillCurrent = () => executionGeneration === networkGeneration && $('network').value === networkId && $('analysis-type').value === 'traffic_growth';


    const button = $('run');

    button.disabled = true;
    button.textContent = 'Analyzing traffic…';
    $('status').textContent = 'Running';

    selectedLinkId = null;
    lastLinks = [];

    $('graph').replaceChildren();
    $('empty').hidden = false;
    $('empty').textContent = 'Evaluating modified traffic demands…';

    $('selected').textContent = 'No link selected.';
    $('scenario-analysis').hidden = true;
    $('scenario-analysis').replaceChildren();

    $('large-network-analysis').hidden = true;
    $('bottleneck-list').replaceChildren();

    // Traffic reports are not ordinary ScenarioRepository reports.
    activeScenarioId = null;
    $('jeylan-question').disabled = true;
    $('jeylan-ask').disabled = true;
    $('ask-selected-link').disabled = true;

    const analysisId = 'traffic-' + crypto.randomUUID().replaceAll('-', '');

    const overrides = [...trafficOverrides].map(([key, factor]) => {
        const [source, target] = key.split(':').map(Number);

        return {source, target, factor};
    });

    try {
        const base = `/api/v1/networks/${encodeURIComponent(networkId)}` + `/analyses/traffic`;

        const workspace = currentWorkspaceSnapshot();

        const response = await fetch(base, {
            method: 'POST', headers: {'Content-Type': 'application/json'}, body: JSON.stringify({
                analysis_id: analysisId,
                demand_model: $('demand').value,
                seed: 42,
                global_scale: workspace.demand_scale,
                commodity_overrides: workspace.commodity_overrides,
                capacity_upgrades: workspace.capacity_upgrades.map(upgrade => ({
                    edge_id: upgrade.edge_id, capacity_factor: upgrade.factor
                }))
            })
        });

        if (!response.ok) {
            const error = await response.json();
            throw new Error(typeof error.detail === 'string' ? error.detail : `Traffic analysis failed (HTTP ${response.status})`);
        }

        const created = await response.json();

        if (created.status !== 'completed') {
            throw new Error('Traffic analysis did not complete.');
        }

        const reportResponse = await fetch(`${base}/${encodeURIComponent(analysisId)}`);

        if (!reportResponse.ok) {
            throw new Error('Unable to retrieve Traffic Growth report.');
        }

        const report = await reportResponse.json();

        if (!stillCurrent()) return;

        renderTrafficReport(report);
        activateJeylanContext(networkId, analysisId, 'traffic_growth');

        $('evidence').textContent = `Traffic analysis: ${analysisId} · ` + `Network: ${networkId} · Mode: fixed_routing`;

    } catch (error) {
        if (stillCurrent()) {
            $('status').textContent = 'Error';
            showError(error.message || String(error));
        }
    } finally {
        if (stillCurrent()) {
            button.disabled = false;
            button.textContent = 'Run traffic analysis';
        }
    }
}