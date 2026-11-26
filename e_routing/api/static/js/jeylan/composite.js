function compositeWorkspaceBrief(workspace) {
    const lines = [];
    const scale = workspace?.demand_scale;
    lines.push(`Demand: ×${Number.isFinite(scale) ? scale : 'N/A'}`);

    const upgrades = Array.isArray(workspace?.capacity_upgrades) ? workspace.capacity_upgrades : [];
    lines.push(upgrades.length ? `Upgrades: ${upgrades.map(item => `edge ${item.edge_id} ×${item.factor}`).join(', ')}` : 'Capacity upgrades: none');

    const overrides = Array.isArray(workspace?.commodity_overrides) ? workspace.commodity_overrides : [];
    lines.push(overrides.length ? `Commodity overrides: ${overrides.length}` : 'Commodity overrides: none');

    return lines;
}

function compositeNumber(value, digits = 2) {
    return typeof value === 'number' && Number.isFinite(value) ? value.toFixed(digits) : 'N/A';
}

function compositePercent(value) {
    if (typeof value !== 'number' || !Number.isFinite(value)) {
        return 'N/A';
    }

    const percent = value * 100;
    return `${percent > 0 ? '+' : ''}${percent.toFixed(2)}%`;
}


function compositeClassForRelativeChange(value) {
    if (typeof value !== 'number' || !Number.isFinite(value)) {
        return 'composite-neutral';
    }
    if (value < 0) return 'composite-good';
    if (value > 0) return 'composite-bad';
    return 'composite-neutral';
}

function compositeFormatValue(value) {
    if (typeof value === 'number' && Number.isFinite(value)) {
        return compositeNumber(value);
    }
    if (typeof value === 'boolean') return value ? 'true' : 'false';
    if (value == null) return 'N/A';
    return typeof value === 'object' ? JSON.stringify(value) : String(value);
}

function compositeElement(tag, text, className = '') {
    const element = document.createElement(tag);
    element.textContent = String(text);

    if (className) {
        element.className = className;
    }

    return element;
}

function clearJeylanCompositeResults() {
    ++compositeTopologyGeneration;
    $('jeylan-composite-results').hidden = true;
    $('composite-status').textContent = '';
    $('composite-summary').textContent = '';
    $('composite-comparisons').replaceChildren();
    $('composite-steps').replaceChildren();
    $('composite-step-detail').replaceChildren();
}

function validateJeylanCompositeResponse(result, networkId) {
    if (result?.network_id !== networkId || result.analysis_type !== 'composite' || result.plan?.kind !== 'composite' || !['completed', 'partial'].includes(result.status) || !Array.isArray(result.completed_steps) || !Array.isArray(result.comparisons) || result.completed_steps.length > 6 || !Array.isArray(result.plan.steps) || result.plan.steps.length < 2 || result.plan.steps.length > 6 || result.workspace_modified !== false) {
        throw new Error('Invalid Jeylan composite response');
    }

    const seenIds = new Set();

    for (const [index, step] of result.completed_steps.entries()) {
        if (step.step_index !== index || step.status !== 'completed' || !COMPOSITE_OPERATION_LABELS[step.operation] || !JEYLAN_ANALYSIS_ID_PATTERN.test(step.analysis_id) || step.operation !== result.plan.steps[index]?.operation || seenIds.has(step.analysis_id) || !step.summary || typeof step.summary !== 'object' || !step.input_workspace || !step.output_workspace) {
            throw new Error('Invalid composite step identity');
        }

        seenIds.add(step.analysis_id);
    }

    if (result.status === 'completed' && result.completed_steps.length !== result.plan.steps.length) {
        throw new Error('Incomplete composite execution');
    }

    for (const comparison of result.comparisons) {
        if (comparison.comparison_type !== 'capacity_n1' || !seenIds.has(comparison.before_analysis_id) || !seenIds.has(comparison.after_analysis_id) || !comparison.counts || typeof comparison.counts !== 'object') {
            throw new Error('Invalid composite comparison identity');
        }
    }
}


function renderCompositeComparisons(result) {
    const root = $('composite-comparisons');
    root.replaceChildren();

    if (result.comparison_status !== 'completed') {
        root.appendChild(compositeElement('p', 'Cross-analysis comparison is not verified or unavailable.'));
        return;
    }

    if (!result.comparisons.length) {
        root.appendChild(compositeElement('p', 'No matching before/after resilience comparison was generated.'));
        return;
    }

    for (const comparison of result.comparisons) {
        const card = compositeElement('div', '', 'composite-comparison');

        card.appendChild(compositeElement('h3', 'Verified N−1 capacity comparison'));

        const metrics = compositeElement('div', '', 'composite-metric-grid');

        const values = [['Before upgrade', compositeNumber(comparison.before_baseline_max_congestion)], ['After upgrade', compositeNumber(comparison.after_baseline_max_congestion)], ['Relative change', compositePercent(comparison.relative_baseline_congestion_change)]];

        for (const [label, value] of values) {
            const metric = compositeElement('div', '', 'composite-metric');

            metric.appendChild(compositeElement('span', label));
            const strong = compositeElement('strong', value);
            strong.className = label === 'Relative change' ? compositeClassForRelativeChange(comparison.relative_baseline_congestion_change) : 'composite-neutral';
            metric.appendChild(strong);
            metrics.appendChild(metric);
        }

        card.appendChild(metrics);

        const counts = comparison.counts;
        const summaryRow = compositeElement('div', '', 'composite-summary-row');

        const chips = [[`${counts.improved} congestion improvements`, counts.improved > 0 ? 'good' : ''], [`${counts.regressed} regressions`, counts.regressed > 0 ? 'bad' : ''], [`${counts.unchanged} unchanged`, ''], [`${counts.became_connected} became connected`, counts.became_connected > 0 ? 'good' : ''], [`${counts.became_disconnected} became disconnected`, counts.became_disconnected > 0 ? 'bad' : ''], [`${counts.disconnected_both} disconnected in both`, '']];

        for (const [label, category] of chips) {
            summaryRow.appendChild(compositeElement('span', label, `composite-chip${category ? ' ' + category : ''}`));
        }

        card.appendChild(summaryRow);
        card.appendChild(compositeElement('div', `Before: ${comparison.before_analysis_id}\n` + `After: ${comparison.after_analysis_id}`, 'composite-comparison-meta'));

        root.appendChild(card);
    }
}


function renderJeylanCompositeResult(result) {
    const root = $('jeylan-composite-results');
    const timeline = $('composite-steps');
    const detail = $('composite-step-detail');

    clearJeylanCompositeResults();

    $('composite-status').textContent = result.status === 'completed' ? 'Completed' : 'Partial';

    $('composite-summary').textContent = `${result.completed_steps.length} of ` + `${result.plan.steps.length} steps completed · ` + `${result.plan.demand_model} demand · ` + 'Simulated workspace only';

    renderCompositeComparisons(result);

    const buttons = [];
    // Each selection has its own token so late report responses cannot
    // overwrite a newer selection or a reset/network switch.
    let selectionToken = 0;
    const compositeNetworkId = result.network_id;
    const compositeNetworkGeneration = networkGeneration;
    const currentCompositeGeneration = compositeTopologyGeneration;

    async function loadStepMeasurements(step, index, token) {
        const current = () => token === selectionToken && currentCompositeGeneration === compositeTopologyGeneration && compositeNetworkGeneration === networkGeneration && $('network').value === compositeNetworkId && !$('jeylan-composite-results').hidden;
        const detail = $('composite-step-detail');
        try {
            const report = await fetchJeylanAnalysisReport(compositeNetworkId, step.operation, step.analysis_id);
            if (!current()) return;
            // This is a view selection, not a workspace change or new run.
            $('analysis-type').value = step.operation;
            await displayJeylanAnalysisResult({
                networkId: compositeNetworkId,
                analysisId: step.analysis_id,
                analysisType: step.operation,
                report,
                stillCurrent: current
            });
            if (!current()) return;
            $('topology-search-feedback').textContent = `Showing persisted step ${index + 1}: ${COMPOSITE_OPERATION_LABELS[step.operation]}.`;
            updateTopologyInspector();
        } catch (error) {
            if (!current()) return;
            $('topology-search-feedback').textContent = `Could not load step ${index + 1} topology: ${error.message}`;
            // Keep the already-persisted timeline and explanatory details.
        }
    }

    function selectStep(index) {
        const step = result.completed_steps[index];
        if (!step) return;

        for (const [i, button] of buttons.entries()) {
            const active = i === index;
            button.classList.toggle('active', active);
            button.setAttribute('aria-pressed', String(active));
        }

        detail.replaceChildren();
        const card = compositeElement('div', '', 'composite-detail-card');
        const header = compositeElement('div', '', 'composite-detail-header');
        header.appendChild(compositeElement('h4', COMPOSITE_OPERATION_LABELS[step.operation]));
        header.appendChild(compositeElement('div', `Analysis ID: ${step.analysis_id}`, 'composite-detail-id'));
        card.appendChild(header);

        const summarySection = compositeElement('section', '', 'composite-detail-section');
        summarySection.appendChild(compositeElement('h5', 'Measured summary'));
        const summaryGrid = compositeElement('div', '', 'composite-kv-grid');
        for (const [key, value] of Object.entries(step.summary)) {
            const item = compositeElement('div', '', 'composite-kv-item');
            item.appendChild(compositeElement('span', key.replaceAll('_', ' ')));
            item.appendChild(compositeElement('strong', compositeFormatValue(value)));
            summaryGrid.appendChild(item);
        }
        summarySection.appendChild(summaryGrid);
        card.appendChild(summarySection);

        const workspaceSection = compositeElement('section', '', 'composite-detail-section');
        workspaceSection.appendChild(compositeElement('h5', 'Workspace transition'));
        const workspaceGrid = compositeElement('div', '', 'composite-workspace-grid');
        for (const [label, workspace] of [['Input workspace', step.input_workspace], ['Output workspace', step.output_workspace]]) {
            const box = compositeElement('div', '', 'composite-workspace-box');
            box.appendChild(compositeElement('h6', label));
            const brief = compositeElement('div', '', 'workspace-brief');
            for (const line of compositeWorkspaceBrief(workspace)) {
                brief.appendChild(compositeElement('div', line));
            }
            box.appendChild(brief);

            const raw = document.createElement('details');
            raw.appendChild(compositeElement('summary', 'View full configuration'));
            const pre = document.createElement('pre');
            pre.textContent = JSON.stringify(workspace, null, 2);
            raw.appendChild(pre);
            box.appendChild(raw);
            workspaceGrid.appendChild(box);
        }
        workspaceSection.appendChild(workspaceGrid);
        card.appendChild(workspaceSection);
        detail.appendChild(card);
        const token = ++selectionToken;
        void loadStepMeasurements(step, index, token);
    }

    for (const [index, step] of result.completed_steps.entries()) {
        const button = compositeElement('button', '', 'composite-step');

        button.type = 'button';

        const content = document.createElement('div');
        content.appendChild(compositeElement('span', `${index + 1}. ` + COMPOSITE_OPERATION_LABELS[step.operation], 'composite-step-title'));

        content.appendChild(compositeElement('span', step.analysis_id, 'composite-step-meta'));

        button.appendChild(content);
        button.appendChild(compositeElement('span', '✓', 'composite-step-status'));

        button.addEventListener('click', () => selectStep(index));

        buttons.push(button);
        timeline.appendChild(button);
    }

    if (result.status === 'partial') {
        const failedIndex = result.failed_step_index;

        timeline.appendChild(compositeElement('p', `Execution stopped at step ${Number.isInteger(failedIndex) ? failedIndex + 1 : 'unknown'}. ${result.error || 'Step failed.'}`, 'composite-bad'));
    }

    root.hidden = false;
    if (buttons.length) {
        selectStep(buttons.length - 1);
    } else {
        detail.textContent = 'No analyses were completed before execution stopped.';
    }
}