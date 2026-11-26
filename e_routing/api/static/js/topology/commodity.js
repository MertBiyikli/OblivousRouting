/** Commodity endpoint overlay. */
(function () {
    let pair = null;
    const $ = id => document.getElementById(id);

    function refreshPair() {
        const svg = $('graph');
        if (!svg) return;
        svg.querySelectorAll('.stitch-commodity-overlay').forEach(e => e.remove());
        if (!pair) return;
        const circles = [...svg.querySelectorAll('circle')];
        const vertices = new Map();
        circles.forEach(c => {
            const t = c.querySelector('title')?.textContent || '';
            const match = /^Vertex\s+(\d+)/.exec(t);
            if (match) vertices.set(Number(match[1]), c);
        });
        const missing = [pair.source, pair.target].filter(id => !vertices.has(id));
        if (missing.length) {
            $('flow-pair-status').textContent = 'Vertex not present in displayed topology: ' + missing.join(', ');
            return;
        }
        const ns = 'http://www.w3.org/2000/svg';
        [[pair.source, 'SOURCE'], [pair.target, 'TARGET']].forEach(([id, role]) => {
            const c = vertices.get(id), x = Number(c.getAttribute('cx')), y = Number(c.getAttribute('cy'));
            const ring = document.createElementNS(ns, 'circle');
            ring.setAttribute('cx', x);
            ring.setAttribute('cy', y);
            ring.setAttribute('r', 12);
            ring.setAttribute('fill', 'none');
            ring.setAttribute('stroke', '#10b981');
            ring.setAttribute('stroke-width', '2');
            ring.setAttribute('class', 'stitch-commodity-overlay');
            ring.style.pointerEvents = 'none';
            svg.append(ring);
            const text = document.createElementNS(ns, 'text');
            text.setAttribute('x', x + 14);
            text.setAttribute('y', y + 4);
            text.setAttribute('class', 'stitch-commodity-overlay stitch-commodity-endpoint-label');
            text.textContent = role + ' ' + id;
            text.style.pointerEvents = 'none';
            svg.append(text);
        });
        $('flow-pair-status').textContent = `Commodity ${pair.source} → ${pair.target}: endpoints highlighted. Routed edge-flow data is not available in this view.`;
    }

    window.stitchRefreshCommodityOverlay = refreshPair;

    document.addEventListener('DOMContentLoaded', () => {
        $('flow-visualize')?.addEventListener('click', () => {
            const a = $('flow-source').value.trim(), b = $('flow-target').value.trim();
            if (!/^\d+$/.test(a) || !/^\d+$/.test(b)) {
                $('flow-pair-status').textContent = 'Enter valid nonnegative source and target vertex IDs.';
                return;
            }
            pair = {source: Number(a), target: Number(b)};
            refreshPair();
        });
        $('flow-clear')?.addEventListener('click', () => {
            pair = null;
            refreshPair();
            $('flow-pair-status').textContent = 'Commodity overlay cleared.';
        });
    });
})();