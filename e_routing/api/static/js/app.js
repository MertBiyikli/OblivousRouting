/**
 * JupITer application bootstrap.
 * Classic script: load all extracted function-declaration files first.
 * Shared state and startup order intentionally preserved.
 */
const $ = id => document.getElementById(id);

$('jeylan-mode-ask').addEventListener('click', () => {
    setJeylanMode('ask');
});

$('jeylan-mode-experiment').addEventListener(
    'click',
    () => setJeylanMode('experiment')
);

const upgrades = new Map();
let lastLinks = [];

// Original physical topology, indexed by network ID.
// Only populate this from a validated baseline report.
const baselineTopologyCache = new Map();

let activeScenarioId = null;
let activeAnalysisContext = null;
let selectedLinkId = null;
let jeylanBusy = false;
let networkGeneration = 0;

// F3-I Topology: one shared graph, one Jeylan instance, no duplicate API state.
const topologySelectedIds = new Set();
let topologyDisplayedLinks = [];
let topologyActive = false;
let compositeTopologyGeneration = 0;
let topologyFailureEdgeId = null;

window.stitchGetTopologyLinks = () => topologyDisplayedLinks.map(link => ({...link}));
window.stitchGetSelectedIds = () => [...topologySelectedIds];

$('nav-topology').addEventListener('click', () => switchStudioView('topology'));
$('nav-analysis').addEventListener('click', () => switchStudioView('analysis'));
$('topology-open-analysis').addEventListener('click', () => switchStudioView('analysis'));

$('topology-inspector-toggle').addEventListener('click', () => {
    const minimized = document.body.classList.toggle('inspector-minimized');
    $('topology-inspector-toggle').textContent = minimized ? '+' : '−';
    $('topology-inspector-toggle').setAttribute('aria-expanded', String(!minimized));
});

$('jeylan-launcher').addEventListener('click', () => {
    setJeylanOpen(true, {focus: true});
});

$('topology-jeylan-toggle').addEventListener('click', () => {
    setJeylanOpen(false, {focus: true});
});

document.addEventListener('keydown', event => {
    if (event.key === 'Escape' &&
        document.body.classList.contains('jeylan-open')) {
        setJeylanOpen(false, {focus: true});
    }
});

setJeylanOpen(false);

window.stitchGetTopologyLinks = () => topologyDisplayedLinks.map(link => ({...link}));
window.stitchGetSelectedIds = () => [...topologySelectedIds];

$('nav-topology').addEventListener(
    'click',
    () => switchStudioView('topology')
);

$('nav-analysis').addEventListener(
    'click',
    () => switchStudioView('analysis')
);

$('topology-clear-selection').addEventListener('click', () => {
    topologySelectedIds.clear();
    selectedLinkId = null;
    $('ask-selected-link').disabled = true;
    syncTopologyEdgeSelection();
});

$('topology-ask-jeylan').addEventListener('click', () => {
    if (!topologySelectedIds.size) {
        $('topology-search-feedback').textContent = 'Select at least one physical edge first.';
        return;
    }
    setJeylanOpen(true);
    setJeylanMode('ask');
    const ids = [...topologySelectedIds].sort((a, b) => a - b);
    $('jeylan-question').value = ids.length === 1
        ? `Analyze physical edge ${ids[0]} using measurements from the displayed persisted analysis. Distinguish load/utilization from verified N-1 failure effects.`
        : `Compare physical edges ${ids.join(', ')} using verified measurements from the displayed persisted analysis. Distinguish load/utilization from measured N-1 failure effects.`;
    if (activeAnalysisContext) {
        $('jeylan-question').focus();
    } else {
        $('topology-search-feedback').textContent = 'An analysis is required before Jeylan can interpret selected edges. You can use Experiment mode to run one.';
    }
});
$('topology-edge-find').addEventListener('click', searchTopologyEdge);
$('topology-edge-search').addEventListener('keydown', event => {
    if (event.key === 'Enter') {
        event.preventDefault();
        searchTopologyEdge();
    }
});


const trafficOverrides = new Map();

$('traffic-add').addEventListener('click', () => {
    const source = Number($('traffic-source').value);
    const target = Number($('traffic-target').value);
    const factor = Number($('traffic-factor').value);

    if (
        $('traffic-source').value.trim() === '' ||
        $('traffic-target').value.trim() === '' ||
        !Number.isSafeInteger(source) ||
        !Number.isSafeInteger(target) ||
        source < 0 ||
        target < 0 ||
        source === target ||
        !Number.isFinite(factor) ||
        factor <= 0 ||
        factor > 1000
    ) {
        showError('Enter distinct nonnegative node IDs and a valid multiplier.');
        return;
    }

    clearError();

    trafficOverrides.set(`${source}:${target}`, factor);

    renderTrafficOverrides();
    refreshWorkspaceSummary();


});


$('jeylan-ask').addEventListener('click', async () => {
    if (jeylanBusy || !activeAnalysisContext) return;

    const question = $('jeylan-question').value.trim();
    if (!question) return;

    const context = activeAnalysisContext;
    const generation = networkGeneration;
    const selectedEdgeIds = topologySelectedIds.size
        ? [...topologySelectedIds].sort((a, b) => a - b)
        : (selectedLinkId === null ? [] : [selectedLinkId]);
    const selectedEdgeId = selectedEdgeIds.length === 1
        ? selectedEdgeIds[0] : null;
    // For N−1, distinguish the currently displayed failure
    // from the physical links selected for comparison.
    const failedEdgeId = context.analysisType === 'n1_resilience'
        ? topologyFailureEdgeId : null;

    const stillCurrent = () =>
        generation === networkGeneration &&
        activeAnalysisContext === context &&
        $('network').value === context.networkId &&
        $('analysis-type').value === context.analysisType;

    jeylanBusy = true;
    $('jeylan-ask').disabled = true;
    $('jeylan-status').textContent =
        'Jeylan is examining the analysis evidence…';

    try {
        const result = await askAnalysisJeylan(
            context,
            question,
            selectedEdgeId,
            selectedEdgeIds,
            failedEdgeId
        );

        if (!stillCurrent()) return;

        renderJeylanAnswer(result.answer);

        const ids = Array.isArray(result.evidence)
            ? result.evidence.map(item => item.id).join(', ')
            : 'Unavailable';

        $('jeylan-status').textContent =
            `Evidence: ${ids} · Analysis: ${context.analysisId}`;

    } catch (error) {
        if (stillCurrent()) {
            $('jeylan-status').textContent =
                error.message || 'Jeylan is unavailable.';
        }

    } finally {
        jeylanBusy = false;

        if (stillCurrent()) {
            $('jeylan-ask').disabled = false;
        }
    }
});


$('jeylan-question').addEventListener('keydown', event => {
    if (event.key === 'Enter' && (event.ctrlKey || event.metaKey)) {
        event.preventDefault();
        $('jeylan-ask').click();
    }
});

$('ask-selected-link').addEventListener('click', () => {
    const link = lastLinks.find(item => item.edge_id === selectedLinkId);
    if (!link || !activeAnalysisContext) return;
    $('jeylan-question').value =
        `Explain this link's recorded load/capacity ratio, whether it is overloaded, ` +
        `what the scenario does and does not establish, and what to investigate next.`;
    $('jeylan-question').focus();
});


$('network').addEventListener('change', resetNetworkState);


// Invalidate previously computed results when the demand model changes.
// Preserve the selected network and all workspace modifications.
$('demand').addEventListener('change', () => {
    // Invalidate in-flight requests for the previous demand model.
    networkGeneration++;
    clearDisplayedResultProvenance();
    clearJeylanCompositeResults();
    // Remove previous analysis results.
    $('scenario-analysis').hidden = true;
    $('scenario-analysis').replaceChildren();

    activeScenarioId = null;
    selectedLinkId = null;
    lastLinks = [];

    // Clear the previous topology visualization.
    $('graph').replaceChildren();
    $('empty').hidden = false;
    $('empty').textContent =
        'Demand model changed. Run an analysis to update the topology.';

    $('large-network-analysis').hidden = true;
    $('bottleneck-list').replaceChildren();
    $('failure-legend').hidden = true;

    $('topology').textContent = 'Network topology';
    $('selected').textContent = 'No link selected.';

    // Clear numerical results from the previous demand model.
    for (const id of [
        'congestion',
        'overloaded',
        'count',
        'baseline-congestion',
        'scenario-congestion',
        'congestion-delta',
        'reference-congestion',
        'upgraded-congestion',
        'capacity-improvement'
    ]) {
        $(id).textContent = '—';
    }

    $('congestion-delta').style.color = '';
    $('capacity-improvement').style.color = '';

    // Previous scenario evidence is no longer current.
    $('evidence').textContent =
        `Demand model changed to ${$('demand').value}. ` +
        'Run an analysis to generate new results.';

    $('interpretation').textContent =
        'The previous results were computed using a different ' +
        'demand model. Run the analysis again.';

    $('status').textContent = 'Ready';

    // Disable Jeylan until new scenario evidence is available.
    $('jeylan-question').value = '';
    $('jeylan-question').disabled = true;
    $('jeylan-ask').disabled = true;
    $('ask-selected-link').disabled = true;
    $('jeylan-status').textContent =
        'Run an analysis before asking Jeylan.';

    clearError();
});


$('add').addEventListener('click', () => {
    const id = Number($('edge').value);

    const factor = Number(
        $('factor').value.trim().replace(',', '.')
    );
    if ($('edge').value === '' || !Number.isInteger(id) || id < 0 || !Number.isFinite(factor) || factor < 1 || factor > 100) {
        showError('Enter a nonnegative integer edge ID and a capacity factor between 1 and 100.');
        return
    }
    clearError();
    upgrades.set(id, factor);
    $('upgrades').textContent = upgrades.size
        ? 'Applied upgrades: ' +
        [...upgrades].map(([edge, factor]) =>
            `Edge ${edge} ×${factor}`).join(' · ')
        : 'Applied upgrades: none';

    refreshWorkspaceSummary();
});


// F3-H: always render the complete supplied topology. Never drop edges
// just because the graph is large; rank lists and drawing are independent.
const LARGE_NETWORK_THRESHOLD = 40;
const INITIAL_RANKING_SIZE = 4;
const RANKING_PAGE_SIZE = 20;
const GRAPH_BOUNDS = {width: 720, height: 340};
let layoutCache = {key: null, positions: null};
let graphView = {x: 0, y: 0, width: 720, height: 340};

$('topology-zoom-in').addEventListener('click', () => zoomGraph(.75,
    graphView.x + graphView.width / 2, graphView.y + graphView.height / 2));
$('topology-zoom-out').addEventListener('click', () => zoomGraph(1.3,
    graphView.x + graphView.width / 2, graphView.y + graphView.height / 2));
$('topology-reset-view').addEventListener('click', resetGraphView);

const graphSvg = $('graph');
const baselineSvg = $('stitch-baseline-graph');

let graphPan = null;

attachGraphNavigation(graphSvg);
attachGraphNavigation(baselineSvg);


$('run').addEventListener('click', async () => {
    let workspace;

    if ($('run').disabled) {
        return;
    }
    try {
        workspace = currentWorkspaceSnapshot();
        validateAnalysisCompatibility(
            $('analysis-type').value,
            workspace
        );
    } catch (error) {
        showError(error.message);
        return;
    }

    // Keep legacy scenario endpoints synchronized with
    // the shared workspace traffic multiplier.
    $('scale').value = String(workspace.demand_scale);
    if ($('analysis-type').value === 'n1_resilience') {
        await executeN1UI();
        return;
    }

    if ($('analysis-type').value === 'link_failure') {
        await executeLinkFailureUI();
        return;
    }

    if ($('analysis-type').value === 'traffic_growth') {
        await executeTrafficGrowthUI();
        return;
    }

    clearError();
    clearDisplayedResultProvenance();

// Capture the validated workspace exactly once.
    const workspaceSnapshot = workspace;

    const scenarioId = 'capacity-' +
        (globalThis.crypto?.randomUUID?.() ||
            String(Date.now()) + Math.random().toString(36).slice(2))
            .replace(/[^A-Za-z0-9_-]/g, '');

// Workspace modifications define the demand baseline.
// Selected capacity upgrades are separate analysis events.
    const payload = {
        scenario_id: scenarioId,
        scenario_kind: 'capacity_upgrade',
        demand_model: $('demand').value,

        workspace: {
            demand_scale: workspaceSnapshot.demand_scale,
            commodity_overrides: workspaceSnapshot.commodity_overrides,
            capacity_upgrades: []
        },

        capacity_upgrades: workspaceSnapshot.capacity_upgrades.map(
            upgrade => ({
                edge_id: upgrade.edge_id,
                capacity_factor: upgrade.factor
            })
        )
    };
    $('run').disabled = true;

    $('run').textContent = 'Simulating…';
    $('status').textContent = 'Running';
    activeScenarioId = null;
    selectedLinkId = null;
    $('ask-selected-link').disabled = true;
    $('selected').textContent = 'Select a link after running a scenario.';
    $('jeylan-question').disabled = true;
    $('jeylan-ask').disabled = true;
    $('jeylan-status').textContent = 'Waiting for scenario execution…';
    $('interpretation').textContent = 'Running deterministic scenario analysis…';

    const executionGeneration = networkGeneration;
    const executionNetworkId = $('network').value;

    if (!executionNetworkId || $('network').disabled) {
        showError('Select a registered network.');
        return;
    }

    const stillCurrent = () =>
        executionGeneration === networkGeneration &&
        $('network').value === executionNetworkId;

    try {
        const executionGeneration = networkGeneration;

// Use the existing selected-network state from Scenario Studio.
// Capture it once so changing the dropdown cannot silently switch networks.
        const networkId = $('network').value;

        if (!networkId) {
            throw new Error("Select a network before running a scenario.");
        }

        const response = await fetch(
            `/api/v1/networks/${encodeURIComponent(networkId)}/scenarios/detailed-analysis`,
            {
                method: 'POST',
                headers: {'Content-Type': 'application/json'},
                body: JSON.stringify(payload)
            }
        );

        const body = await response.json();

        if (!response.ok) {
            throw new Error(
                typeof body.detail === 'string'
                    ? body.detail
                    : JSON.stringify(body.detail || body)
            );
        }

        if (!stillCurrent()) return;


        if (
            body.status !== 'completed' ||
            !body.report_available ||
            body.analysis?.scenario_kind !== 'capacity_upgrade'
        ) {
            throw new Error('Incomplete capacity analysis response');
        }

        await displayCapacityAnalysisResult({
            networkId,
            scenarioId: body.scenario_id,
            report: body.analysis,
            stillCurrent
        });

    } catch (e) {
        if (stillCurrent()) {
            $('status').textContent = 'Error';
            showError(e.message || String(e));
        }
    } finally {
        if (stillCurrent()) {
            $('run').disabled = false;
            $('run').textContent =
                $('analysis-type').value === 'link_failure'
                    ? 'Run failure analysis'
                    : 'Run simulation';
        }
    }

});


// initial load
$('run').disabled = true;
loadNetworks();

$('reset-network').addEventListener(
    'click', resetWorkspace
);


$('traffic-scale').addEventListener(
    'input', refreshWorkspaceSummary
);

$('analysis-type').addEventListener(
    'change', refreshWorkspaceSummary
);

$('analysis-type').addEventListener(
    'change',
    updateAnalysisMode
);

$('import-network').addEventListener('click', async () => {
    const name = $('import-name').value.trim();
    const file = $('import-file').files[0];

    if (!name || !file) {
        showError('Enter a network name and select an LGF file.');
        return;
    }

    if (!file.name.toLowerCase().endsWith('.lgf')) {
        showError('Only LGF graph files are supported.');
        return;
    }

    if (file.size > 8 * 1024 * 1024) {
        showError('Graph file exceeds the 8 MiB upload limit.');
        return;
    }

    const button = $('import-network');
    const status = $('import-status');

    button.disabled = true;
    button.textContent = 'Importing network…';
    status.textContent = 'Uploading graph and generating routing scheme…';
    clearError();

    try {
        const data = new FormData();
        data.append('name', name);
        data.append('file', file);

        const response = await fetch('/api/v1/networks/import', {
            method: 'POST',
            body: data
        });

        const result = await response.json();

        if (!response.ok) {
            throw new Error(
                typeof result.detail === 'string'
                    ? result.detail
                    : `Import failed (HTTP ${response.status})`
            );
        }

        if (result.status !== 'ready' || !result.network_id) {
            throw new Error('Invalid network import response.');
        }

        await loadNetworks(result.network_id);

        if ($('network').value !== result.network_id) {
            throw new Error(
                'Imported network could not be activated.'
            );
        }

        $('network').dispatchEvent(new Event('change'));

        $('network-info').textContent =
            `Active network: ${result.name} · ` +
            `${result.nodes} nodes · ${result.edges} edges`;

        $('status').textContent = 'Ready';
        $('evidence').textContent =
            `Network ${result.name} imported and selected. ` +
            'Run an analysis to generate results.';

        $('interpretation').textContent =
            `The active network is now ${result.name}. ` +
            'Configure and execute an analysis.';

        status.textContent =
            `Imported ${result.name}: ` +
            `${result.nodes} nodes, ${result.edges} edges.`;

        $('import-file').value = '';

    } catch (error) {
        status.textContent = 'Network import failed.';
        showError(error.message || String(error));

    } finally {
        button.disabled = false;
        button.textContent = 'Import network';
    }
});


const JEYLAN_ANALYSIS_ID_PATTERN = /^jeylan-[a-f0-9]{32}$/;

const COMPOSITE_OPERATION_LABELS = {
    traffic_growth: 'Traffic Growth',
    capacity_upgrade: 'Capacity Upgrade',
    n1_resilience: 'N−1 Resilience',
    link_failure: 'Link Failure'
};


const JEYLAN_ANALYSIS_TYPES = new Set([
    'capacity_upgrade',
    'link_failure',
    'n1_resilience',
    'traffic_growth'
]);


$('jeylan-run-experiment').addEventListener('click', async () => {
    if (jeylanBusy) return;

    const networkId = $('network').value;
    const question =
        $('jeylan-experiment-question').value.trim();

    if (!networkId || $('network').disabled) {
        $('jeylan-experiment-status').textContent =
            'Select a registered network first.';
        return;
    }

    if (!question) {
        $('jeylan-experiment-status').textContent =
            'Describe the analysis to execute.';
        return;
    }

    let workspace;
    try {
        workspace = currentWorkspaceSnapshot();
    } catch (error) {
        $('jeylan-experiment-status').textContent =
            error.message;
        return;
    }

    const generation = networkGeneration;
    const selectedMode = $('analysis-type').value;
    const originalWorkspace = JSON.stringify(workspace);

    const stillCurrent = () =>
        generation === networkGeneration &&
        $('network').value === networkId &&
        $('analysis-type').value === selectedMode;

    jeylanBusy = true;
    $('jeylan-run-experiment').disabled = true;
    $('jeylan-experiment-status').textContent =
        'Jeylan is executing a simulated network analysis…';

    try {
        const response = await fetch(
            `/api/v1/networks/${encodeURIComponent(networkId)}` +
            '/assistant/analyses/ask',
            {
                method: 'POST',
                headers: {'Content-Type': 'application/json'},
                body: JSON.stringify({question, workspace})
            }
        );

        const result = await response.json();

        if (!response.ok) {
            throw new Error(
                typeof result.detail === 'string'
                    ? result.detail
                    : JSON.stringify(
                        result.detail ||
                        `Analysis failed (HTTP ${response.status})`
                    )
            );
        }


        if (!stillCurrent()) return;

        if (result.analysis_type === 'composite') {
            validateJeylanCompositeResponse(result, networkId);

            if (
                JSON.stringify(currentWorkspaceSnapshot()) !==
                originalWorkspace
            ) {
                throw new Error(
                    'Workspace changed during Jeylan composite execution'
                );
            }

            renderJeylanCompositeResult(result);
            renderJeylanAnswer(result.answer, result);

            $('jeylan-experiment-status').textContent =
                `${result.status === 'completed'
                    ? 'Completed'
                    : 'Partially completed'} · ` +
                `${result.completed_steps.length} steps`;

            $('evidence').textContent =
                `Jeylan composite investigation · ` +
                `${result.completed_steps.length} persisted analyses · ` +
                `Workspace unchanged`;

            return;
        }

// Preserve the existing single-analysis path.
        if (
            result.network_id !== networkId ||
            result.status !== 'completed' ||
            !JEYLAN_ANALYSIS_TYPES.has(result.analysis_type) ||
            !JEYLAN_ANALYSIS_ID_PATTERN.test(result.analysis_id) ||
            result.plan?.kind !== result.analysis_type
        ) {
            throw new Error('Invalid Jeylan analysis identity');
        }

        clearJeylanCompositeResults();


        const report = await fetchJeylanAnalysisReport(
            networkId,
            result.analysis_type,
            result.analysis_id
        );

        if (!stillCurrent()) return;

        // Analysis type is a display selection, not a new
        // workspace modification. Do not dispatch a change
        // event that would reset the displayed analysis.
        $('analysis-type').value = result.analysis_type;

        const resultStillCurrent = () =>
            generation === networkGeneration &&
            $('network').value === networkId &&
            $('analysis-type').value === result.analysis_type;

        await displayJeylanAnalysisResult({
            networkId,
            analysisId: result.analysis_id,
            analysisType: result.analysis_type,
            report,
            stillCurrent: resultStillCurrent
        });

        if (!resultStillCurrent()) return;

        // The backend remains the source of truth.
        // The explanation is displayed alongside the report.
        renderJeylanAnswer(result.answer);

        $('jeylan-experiment-status').textContent =
            `Completed · ${result.analysis_type} · ` +
            result.analysis_id;

        // The editor must not silently acquire the experiment
        // changes. Only the displayed analysis changes.
        if (
            JSON.stringify(currentWorkspaceSnapshot()) !==
            originalWorkspace
        ) {
            throw new Error(
                'Workspace changed during Jeylan analysis'
            );
        }

    } catch (error) {
        if (
            generation === networkGeneration &&
            $('network').value === networkId
        ) {
            $('jeylan-experiment-status').textContent =
                error.message || 'Jeylan analysis unavailable.';
        }
    } finally {
        jeylanBusy = false;
        $('jeylan-run-experiment').disabled = false;
    }
});


updateAnalysisMode();