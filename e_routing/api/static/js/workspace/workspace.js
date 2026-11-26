function currentWorkspaceSnapshot() {
    const demandScale = Number($('traffic-scale').value);

    if (!Number.isFinite(demandScale) || demandScale <= 0 || demandScale > 100) {
        throw new Error('Demand multiplier must be greater than 0 and at most 100.');
    }

    const capacity_upgrades = [...upgrades]
        .map(([edge_id, factor]) => ({
            edge_id: Number(edge_id), factor: Number(factor)
        }))
        .sort((a, b) => a.edge_id - b.edge_id);

    const commodity_overrides = [...trafficOverrides]
        .map(([key, factor]) => {
            const [source, target] = key.split(':').map(Number);
            return {source, target, factor: Number(factor)};
        })
        .sort((a, b) => a.source - b.source || a.target - b.target);

    return {
        capacity_upgrades, demand_scale: demandScale, commodity_overrides
    };
}

function workspaceIsModified(workspace) {
    return (workspace.capacity_upgrades.length > 0 || workspace.demand_scale !== 1 || workspace.commodity_overrides.length > 0);
}

function refreshWorkspaceSummary() {
    const summary = $('workspace-summary');

    try {
        const workspace = currentWorkspaceSnapshot();

        if (!workspaceIsModified(workspace)) {
            summary.textContent = 'Original network · No modifications';
            return;
        }

        summary.textContent = [`${workspace.capacity_upgrades.length} capacity upgrades`, `Traffic ×${workspace.demand_scale}`, `${workspace.commodity_overrides.length} commodity overrides`].join(' · ');
    } catch (error) {
        summary.textContent = error.message;
    }
}


function resetWorkspace() {
    clearDisplayedResultProvenance();
    clearJeylanCompositeResults();
    // Invalidate results belonging to an earlier workspace.
    networkGeneration++;

    upgrades.clear();
    trafficOverrides.clear();

    $('traffic-scale').value = '1.0';
    $('scale').value = '1.0';

    $('traffic-source').value = '';
    $('traffic-target').value = '';
    $('traffic-factor').value = '2.0';

    $('edge').value = '';
    $('factor').value = '1.5';
    $('failed-edge').value = '';

    selectedLinkId = null;
    activeScenarioId = null;
    lastLinks = [];

    $('graph').replaceChildren();
    $('scenario-analysis').replaceChildren();
    $('scenario-analysis').hidden = true;
    $('failure-legend').hidden = true;
    $('large-network-analysis').hidden = true;
    $('bottleneck-list').replaceChildren();

    $('topology').textContent = 'Network topology';
    $('status').textContent = 'Ready';

    for (const id of ['congestion', 'overloaded', 'count', 'baseline-congestion', 'scenario-congestion', 'congestion-delta', 'reference-congestion', 'upgraded-congestion', 'capacity-improvement']) {
        $(id).textContent = '—';
    }

    $('congestion-delta').style.color = '';
    $('capacity-improvement').style.color = '';

    $('empty').hidden = false;
    $('empty').textContent = 'Workspace reset. Run an analysis to load the original network.';

    $('evidence').textContent = 'Workspace reset · No analysis executed.';

    $('interpretation').textContent = 'The original workspace has been restored. ' + 'Run an analysis to display new results.';

    $('jeylan-question').value = '';
    $('jeylan-status').textContent = 'Run an analysis before asking Jeylan.';

    $('jeylan-question').disabled = true;
    $('jeylan-ask').disabled = true;
    $('ask-selected-link').disabled = true;

    renderTrafficOverrides();

    // Use the existing capacity-list renderer here if available.
    $('upgrades').textContent = 'Applied upgrades: none';

    clearError();
    refreshWorkspaceSummary();
    invalidateJeylanContext();

    $('run').disabled = !$('network').value;
}

function validateAnalysisCompatibility(type, workspace) {
    if (type === 'n1_resilience') {
        if ($('execution-mode').value === 'research') {
            throw new Error('N−1 RESEARCH mode is not supported by the API yet.');
        }
        return;
    }


    if (type === 'capacity_upgrade' && workspace.capacity_upgrades.length < 1) {
        throw new Error('Capacity Planning requires at least one capacity upgrade. ' + 'Enter an edge ID and capacity factor, then click Add capacity upgrade.');
    }
}