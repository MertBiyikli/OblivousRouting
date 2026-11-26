async function loadNetworks(selectedNetworkId = null) {
    const selector = $('network');
    selector.disabled = true;
    selector.replaceChildren(new Option('Loading networks...', ''));

    try {
        const response = await fetch('/api/v1/networks', {
            cache: 'no-store'
        });

        if (!response.ok) {
            throw new Error(`Network registry returned HTTP ${response.status}`);
        }

        const data = await response.json();

        if (!Array.isArray(data.networks)) {
            throw new Error('Invalid network registry response');
        }

        selector.replaceChildren();

        for (const network of data.networks) {
            selector.add(new Option(`${network.name} · ${network.solver}`, network.id));
        }

        if (data.networks.length === 0) {
            selector.add(new Option('No registered networks', ''));
            $('network-info').textContent = 'No networks are currently registered.';
            $('run').disabled = true;
            return;
        }

        selector.disabled = false;

        if (selectedNetworkId) {
            const exists = data.networks.some(network => network.id === selectedNetworkId);

            if (!exists) {
                throw new Error('Imported network was not found in the registry.');
            }

            selector.value = selectedNetworkId;
            resetNetworkState();
        }

        $('run').disabled = false;

        $('network-info').textContent = `${data.networks.length} registered network(s) available.`;

    } catch (error) {
        selector.replaceChildren(new Option('Network registry unavailable', ''));

        $('run').disabled = true;
        $('network-info').textContent = 'Unable to load registered networks.';

        showError(error.message);
    }
}

function resetNetworkState() {
    clearDisplayedResultProvenance();
    clearJeylanCompositeResults();
    $('scenario-analysis').hidden = true;
    $('scenario-analysis').replaceChildren();
    $('run').textContent = 'Run simulation';
    $('run').disabled = !$('network').value || $('network').disabled;
    // Invalidate results and topology from the previous network.
    upgrades.clear();
    trafficOverrides.clear();
    $('traffic-scale').value = '1.0';
    $('scale').value = '1.0';

    renderTrafficOverrides();
    refreshWorkspaceSummary();
    lastLinks = [];
    activeScenarioId = null;
    selectedLinkId = null;
    networkGeneration++;

    $('edge').value = '';
    $('upgrades').textContent = 'Applied upgrades: none';

    $('graph').replaceChildren();
    $('failure-legend').hidden = true;
    $('selected').textContent = 'Select a link after running a scenario.';

    $('topology').textContent = 'Network topology';
    $('empty').hidden = false;
    $('empty').textContent = 'Run a scenario to load the selected network topology.';

    $('large-network-analysis').hidden = true;
    $('bottleneck-list').replaceChildren();
    $('topology-summary').textContent = '';


    const metricIds = ['congestion', 'overloaded', 'count', 'baseline-congestion', 'scenario-congestion', 'congestion-delta', 'reference-congestion', 'upgraded-congestion'];

    for (const id of metricIds) {
        const element = $(id);
        if (element) element.textContent = '—';
    }

    lastLinks = [];
    activeScenarioId = null;
    selectedLinkId = null;

    $('capacity-improvement').textContent = '—';
    $('evidence').textContent = '';
    $('interpretation').textContent = 'Run a scenario to inspect the selected network.';

    $('jeylan-question').value = '';
    $('jeylan-question').disabled = true;
    $('jeylan-ask').disabled = true;
    $('ask-selected-link').disabled = true;
    $('jeylan-status').textContent = 'Run a scenario before asking Jeylan.';

    $('status').textContent = 'Ready';
    for (const id of ['congestion', 'overloaded', 'count', 'baseline-congestion', 'scenario-congestion', 'congestion-delta', 'reference-congestion', 'upgraded-congestion', 'capacity-improvement']) {
        $(id).textContent = '—';
    }

    $('congestion-delta').style.color = '';
    $('capacity-improvement').style.color = '';

    $('evidence').textContent = 'Workspace reset · No analysis executed.';

    $('interpretation').textContent = 'The original network workspace has been restored. ' + 'Run an analysis to display new results.';

    invalidateJeylanContext();
    clearError();
}