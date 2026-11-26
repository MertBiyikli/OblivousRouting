function topologyLayout(links) {
    const vertices = [...new Set(links.flatMap(l => [l.source, l.target]))]
        .sort((a, b) => Number(a) - Number(b));
    const key = vertices.join(',') + '|' + links
        .map(l => `${Math.min(l.source, l.target)}:${Math.max(l.source, l.target)}`)
        .sort().join(',');
    if (layoutCache.key === key) return layoutCache.positions;
    const positions = new Map();
    if (!vertices.length) return positions;
    const adjacency = new Map(vertices.map(v => [v, new Set()]));
    for (const l of links) {
        adjacency.get(l.source)?.add(l.target);
        adjacency.get(l.target)?.add(l.source);
    }
    // Deterministic graph-distance layering: physically connected neighbors
    // determine placement, not vertex ID order around an ellipse.
    const seen = new Set();
    const components = [];
    for (const vertex of vertices) {
        if (seen.has(vertex)) continue;
        const component = [];
        const queue = [vertex];
        seen.add(vertex);
        for (let i = 0; i < queue.length; i++) {
            const v = queue[i];
            component.push(v);
            for (const neighbor of [...adjacency.get(v)].sort((a, b) => a - b)) {
                if (!seen.has(neighbor)) {
                    seen.add(neighbor);
                    queue.push(neighbor);
                }
            }
        }
        components.push(component);
    }
    components.sort((a, b) => b.length - a.length);
    let yOffset = 0;
    for (const component of components) {
        // A high-degree vertex is a useful deterministic BFS root.
        const root = [...component].sort((a, b) => adjacency.get(b).size - adjacency.get(a).size || a - b)[0];
        const layers = [];
        const depth = new Map([[root, 0]]);
        const queue = [root];
        for (let i = 0; i < queue.length; i++) {
            const v = queue[i], d = depth.get(v);
            (layers[d] ||= []).push(v);
            for (const u of [...adjacency.get(v)].sort((a, b) => a - b)) {
                if (!depth.has(u)) {
                    depth.set(u, d + 1);
                    queue.push(u);
                }
            }
        }
        const maxLayer = Math.max(1, layers.length - 1);
        for (let d = 0; d < layers.length; d++) {
            const layer = layers[d];
            for (let j = 0; j < layer.length; j++) {
                positions.set(layer[j], {
                    x: 48 + 624 * d / maxLayer, y: 26 + yOffset + 265 * (j + 1) / (layer.length + 1)
                });
            }
        }
        yOffset += components.length > 1 ? 290 : 0;
    }
    // Iterative neighbor attraction + short-range repulsion spreads clusters.
    // Cap quadratic repulsion work for very large imported networks.
    if (vertices.length <= 260 && components.length === 1) {
        for (let iteration = 0; iteration < 75; iteration++) {
            const forces = new Map(vertices.map(v => [v, {x: 0, y: 0}]));
            for (const l of links) {
                const a = positions.get(l.source), b = positions.get(l.target);
                if (!a || !b) continue;
                const dx = b.x - a.x, dy = b.y - a.y;
                const distance = Math.max(1, Math.hypot(dx, dy));
                const strength = (distance - 65) * .007;
                const fx = dx / distance * strength, fy = dy / distance * strength;
                forces.get(l.source).x += fx;
                forces.get(l.source).y += fy;
                forces.get(l.target).x -= fx;
                forces.get(l.target).y -= fy;
            }
            for (let i = 0; i < vertices.length; i++) for (let j = i + 1; j < vertices.length; j++) {
                const a = positions.get(vertices[i]), b = positions.get(vertices[j]);
                const dx = a.x - b.x, dy = a.y - b.y, dist2 = Math.max(16, dx * dx + dy * dy);
                if (dist2 > 130 * 130) continue;
                const magnitude = 180 / dist2;
                const fx = dx * magnitude, fy = dy * magnitude;
                forces.get(vertices[i]).x += fx;
                forces.get(vertices[i]).y += fy;
                forces.get(vertices[j]).x -= fx;
                forces.get(vertices[j]).y -= fy;
            }
            for (const v of vertices) {
                const p = positions.get(v), f = forces.get(v);
                p.x = Math.max(18, Math.min(702, p.x + Math.max(-8, Math.min(8, f.x))));
                p.y = Math.max(16, Math.min(324, p.y + Math.max(-8, Math.min(8, f.y))));
            }
        }
    }
    // Fit complete layout, including disconnected components, into the canvas.
    const xs = [...positions.values()].map(p => p.x);
    const ys = [...positions.values()].map(p => p.y);
    const minX = Math.min(...xs), maxX = Math.max(...xs), minY = Math.min(...ys), maxY = Math.max(...ys);
    const scale = Math.min(660 / Math.max(1, maxX - minX), 280 / Math.max(1, maxY - minY));
    const offsetX = (720 - (maxX - minX) * scale) / 2, offsetY = (340 - (maxY - minY) * scale) / 2;
    for (const p of positions.values()) {
        p.x = offsetX + (p.x - minX) * scale;
        p.y = offsetY + (p.y - minY) * scale;
    }
    layoutCache = {key, positions};
    return positions;
}

function selectTopologyLinks(links) {
    const vertices = new Set(links.flatMap(l => [l.source, l.target]));
    return {
        large: vertices.size > LARGE_NETWORK_THRESHOLD || links.length > 150,
        visible: links,
        ranked: [...links].sort((a, b) => (Number(b.utilization) || 0) - (Number(a.utilization) || 0))
    };
}

function drawGraph(links, options = {}) {
    const svg = $('graph');
    svg.replaceChildren();
    resetGraphView();
    // Baseline links preserve physical topology and the position of failed edges.
    // The inspector must use the measured *current* links for survivors.
    const topologyLinks = options.baselineLinks || links;
    window.stitchTopologyBaseline = topologyLinks.map(link => {
        const baselineUtil = Number(link.baseline_utilization);
        return {
            ...link,
            utilization: link.baseline_utilization != null && Number.isFinite(baselineUtil) ? baselineUtil : link.utilization,
            capacity: link.baseline_capacity ?? link.capacity,
            load: link.baseline_load ?? link.load
        };
    });
    window.stitchTopologyCurrent = links.map(link => ({...link}));
    window.stitchTopologyFailed = [...(options.failedIds || [])];
    const measuredById = new Map(links.map(link => [link.edge_id, link]));
    const failedIds = new Set(options.failedIds || []);
    topologyFailureEdgeId = failedIds.size === 1 ? [...failedIds][0] : null;
    topologyDisplayedLinks = topologyLinks.map(link => failedIds.has(link.edge_id) ? {
        ...link,
        topology_failed: true
    } : (measuredById.get(link.edge_id) || {...link, topology_unmeasured: true}));
    // Keep selections only while their links remain present in the displayed graph.
    for (const id of topologySelectedIds) {
        if (!topologyLinks.some(link => link.edge_id === id)) topologySelectedIds.delete(id);
    }
    updateTopologyInspector();
    const selection = selectTopologyLinks(topologyLinks);
    const analysis = $('large-network-analysis');
    analysis.hidden = !selection.large;
    $('topology-mode').textContent = selection.large ? 'Topology + bottleneck ranking' : '';
    if (selection.large) {
        $('topology-summary').textContent = `All ${topologyLinks.length} physical links drawn. ` + `Use zoom and pan to explore connectivity. Ranked results are listed below.`;
    }
    if (!topologyLinks.length) return;
    const positions = topologyLayout(topologyLinks);
    const vertices = [...positions.keys()];
    const failed = topologyLinks.filter(l => failedIds.has(l.edge_id));
    const survivors = links.filter(l => !failedIds.has(l.edge_id));

    function addLine(l, isFailed) {
        const a = positions.get(l.source), b = positions.get(l.target);
        if (!a || !b) return;
        const util = Number(l.utilization);
        const color = isFailed ? '#ffb454' : util > 2 ? '#f26873' : util > 1 ? '#f0aa36' : '#40b7ef';
        const line = svgEl('line', {
            x1: a.x,
            y1: a.y,
            x2: b.x,
            y2: b.y,
            stroke: color,
            'stroke-width': isFailed ? 3 : 2.2,
            'stroke-opacity': isFailed ? .95 : .75, ...(isFailed ? {'stroke-dasharray': '8 6'} : {})
        });
        if (isFailed) line.dataset.failedEdgeId = String(l.edge_id); else line.dataset.edgeId = String(l.edge_id);

        if (l.edge_id === selectedLinkId && !isFailed) line.classList.add('selected-edge');
        line.addEventListener('click', (event) => {
            if (topologyActive) {
                selectTopologyEdge(l, Boolean(event.shiftKey || event.metaKey || event.ctrlKey));
                return;
            }
            topologySelectedIds.clear();
            if (!isFailed) topologySelectedIds.add(l.edge_id);
            updateTopologyInspector();
            svg.querySelectorAll('line.selected-edge').forEach(el => el.classList.remove('selected-edge'));
            selectedLinkId = isFailed ? null : l.edge_id;
            if (!isFailed) line.classList.add('selected-edge');
            $('edge').value = l.edge_id;
            $('failed-edge').value = l.edge_id;
            $('selected').textContent = isFailed ? `FAILED physical link ${l.edge_id} (${l.source} → ${l.target}) · removed from modified network` : `Selected edge ${l.edge_id} (${l.source} → ${l.target}) · capacity ${l.capacity} · utilization ${Number.isFinite(util) ? util.toFixed(3) : 'N/A'}`;
            $('ask-selected-link').disabled = isFailed || !activeAnalysisContext;
            if (!isFailed && activeAnalysisContext) $('jeylan-status').textContent = `Edge ${l.edge_id} selected · Ready to ask Jeylan`;
        });
        svg.appendChild(line);
    }

    for (const l of survivors) addLine(l, false);
    for (const l of failed) addLine(l, true);
    const large = vertices.length > 40;
    for (const v of vertices) {
        const p = positions.get(v);
        const circle = svgEl('circle', {cx: p.x, cy: p.y, r: large ? 3.3 : 6});
        const title = svgEl('title', {});
        title.textContent = `Vertex ${v}`;
        circle.appendChild(title);
        svg.appendChild(circle);
        if (!large) {
            const t = svgEl('text', {x: p.x + 8, y: p.y - 7});
            t.textContent = String(v);
            svg.appendChild(t);
        }
    }
    $('empty').hidden = true;
    syncTopologyEdgeSelection();
    if (window.stitchRefreshCommodityOverlay) window.stitchRefreshCommodityOverlay();
    if (window.stitchRefreshTopology) window.stitchRefreshTopology();
}