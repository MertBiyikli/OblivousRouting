function topologyFormat(value, decimals = 3) {
    return typeof value === 'number' && Number.isFinite(value)
        ? value.toFixed(decimals) : 'N/A';
}

function updateTopologyInspector() {
    const links = topologyDisplayedLinks;
    const n = new Set(links.flatMap(edge => [edge.source, edge.target])).size;
    const overloads = links.filter(edge =>
        !edge.topology_failed && !edge.topology_unmeasured &&
        Number(edge.utilization) > 1
    ).length;
    const summary = $('topology-inspector-summary');
    const baseline = $('baseline-congestion').textContent.trim();
    const scenario = $('scenario-congestion').textContent.trim();
    const composite = $('jeylan-composite-results');
    const compositeSummary = composite && !composite.hidden
        ? $('composite-summary').textContent.trim() : '';

    const detail = $('topology-edge-details');
    detail.replaceChildren();
    const selected = links.filter(edge => topologySelectedIds.has(edge.edge_id));

    const total = document.createElement('div');
    total.className = 'topology-hint';
    total.textContent = selected.length === 1
        ? 'Selected physical link'
        : `${selected.length} selected physical links · ${selected.filter(edge => Number(edge.utilization) > 1).length} overloaded`;
    detail.appendChild(total);
    for (const edge of selected.slice(0, 25)) {
        const card = document.createElement('div');
        card.className = 'topology-edge-card';
        const heading = document.createElement('strong');
        heading.textContent = `Edge ${edge.edge_id} · ${edge.source} → ${edge.target}`;
        const content = document.createElement('div');
        content.textContent = edge.topology_failed
            ? 'Failed physical link · Removed from this scenario; post-failure load and utilization are not defined.'
            : edge.topology_unmeasured
                ? 'No modified measurement available for this physical link.'
                : `Capacity ${topologyFormat(edge.capacity)} · Load ${topologyFormat(edge.load)} · Utilization ${topologyFormat(edge.utilization)}×`;
        card.append(heading, content);
        detail.appendChild(card);
    }
    if (selected.length > 25) {
        const note = document.createElement('div');
        note.className = 'topology-hint';
        note.textContent = `Showing 25 of ${selected.length} selected links.`;
        detail.appendChild(note);
    }
}

function syncTopologyEdgeSelection() {
    for (const line of $('graph').querySelectorAll('line[data-edge-id], line[data-failed-edge-id]')) {
        const id = Number(line.dataset.edgeId ?? line.dataset.failedEdgeId);
        line.classList.toggle('topology-multi-selected', topologySelectedIds.has(id));
    }
    updateTopologyInspector();
}

function selectTopologyEdge(link, additive = false) {
    if (!additive) topologySelectedIds.clear();
    if (additive && topologySelectedIds.has(link.edge_id)) topologySelectedIds.delete(link.edge_id);
    else topologySelectedIds.add(link.edge_id);
    // Jeylan's existing single-edge question continues to use the primary selected ID.
    selectedLinkId = topologySelectedIds.size
        ? [...topologySelectedIds][topologySelectedIds.size - 1] : null;
    if (selectedLinkId !== null) {
        const primary = topologyDisplayedLinks.find(edge => edge.edge_id === selectedLinkId) || link;
        $('edge').value = selectedLinkId;
        $('failed-edge').value = selectedLinkId;
        $('selected').textContent = `Selected edge ${selectedLinkId} (${primary.source} → ${primary.target})`;
    }
    $('ask-selected-link').disabled = !activeAnalysisContext || selectedLinkId === null;
    syncTopologyEdgeSelection();
}

function searchTopologyEdge() {
    const query = $('topology-edge-search').value.trim();
    if (!query) {
        $('topology-search-feedback').textContent = 'Enter a physical edge ID or two endpoints.';
        return;
    }
    const values = (query.match(/\d+/g) || []).map(Number);
    let matches = [];
    if (values.length === 1) {
        matches = topologyDisplayedLinks.filter(edge => edge.edge_id === values[0]);
    } else if (values.length === 2) {
        matches = topologyDisplayedLinks.filter(edge =>
            (edge.source === values[0] && edge.target === values[1]) ||
            (edge.source === values[1] && edge.target === values[0]));
    }
    if (!matches.length) {
        $('topology-search-feedback').textContent = 'No matching physical link in the displayed topology.';
        return;
    }
    topologySelectedIds.clear();
    for (const link of matches) topologySelectedIds.add(link.edge_id);
    selectedLinkId = matches[0].edge_id;
    $('edge').value = selectedLinkId;
    $('failed-edge').value = selectedLinkId;
    $('ask-selected-link').disabled = !activeAnalysisContext;
    $('topology-search-feedback').textContent = `${matches.length} matching link(s) selected. Zoom in to inspect them.`;
    syncTopologyEdgeSelection();
}