/**
 * Topology Workbench
 *
 * Responsibilities:
 * - Single, Compare and Difference visualization modes.
 * - Topology filtering.
 * - Network telemetry.
 * - Edge audit table.
 * - Synchronization with the existing graph renderer.
 *
 * Public interface:
 * - window.stitchRefreshTopology()
 *
 * Existing topology state and graph ownership remain unchanged.
 */

(function () {
    let mode = 'single', filter = 'all';

    const $ = id => document.getElementById(id);
    const num = x => Number.isFinite(Number(x)) ? Number(x) : null;

    function links() {
        return window.stitchGetTopologyLinks?.() || [];
    }

    function paintMode(svg, all) {
        if (!document.body.classList.contains('topology-view')) return;

        document.body.classList.toggle(
            'stitch-mode-compare',
            mode === 'compare'
        );

        document.body.classList.toggle(
            'stitch-mode-difference',
            mode === 'difference'
        );

        const baseline = window.stitchTopologyBaseline || [];

        const baselineById = new Map(
            baseline.map(l => [String(l.edge_id), l])
        );

        const currentById = new Map(
            (window.stitchTopologyCurrent || all)
                .map(l => [String(l.edge_id), l])
        );

        const failed = new Set(
            (window.stitchTopologyFailed || []).map(String)
        );

        const baselineSvg = $('stitch-baseline-graph');

        baselineSvg.replaceChildren();

        baselineSvg.setAttribute(
            'viewBox',
            svg.getAttribute('viewBox') || '0 0 720 340'
        );

        if (mode === 'compare') {
            // Clone exact geometry, labels, and vertices from
            // the live graph; then restyle baseline edges using
            // baseline measurements, not scenario values.

            for (const child of svg.children) {
                if (
                    child.classList?.contains('stitch-commodity-overlay') ||
                    child.classList?.contains('stitch-delta-label')
                ) continue;

                const clone = child.cloneNode(true);

                if (clone.tagName.toLowerCase() === 'line') {
                    const id =
                        clone.dataset.edgeId ||
                        clone.dataset.failedEdgeId;

                    if (!baselineById.has(String(id))) continue;

                    const baselineEdge = baselineById.get(String(id));
                    const util = Number(baselineEdge.utilization);

                    clone.style.stroke =
                        Number.isFinite(util) && util > 1
                            ? '#b14b4b'
                            : '#9a989f';

                    clone.style.strokeWidth =
                        Number.isFinite(util) && util > 1
                            ? '3'
                            : '2';

                    clone.style.strokeDasharray = '';
                    clone.style.opacity = '1';

                    clone.classList.remove('selected-edge');

                    clone.classList.toggle(
                        'stitch-filtered-out',
                        child.classList.contains('stitch-filtered-out')
                    );
                }

                baselineSvg.append(clone);
            }
        }

        svg.querySelectorAll('.stitch-delta-label')
            .forEach(el => el.remove());

        svg.querySelectorAll('line').forEach(line => {
            const classes = [
                'stitch-delta-improved',
                'stitch-delta-worsened',
                'stitch-delta-strong-improvement',
                'stitch-delta-improvement',
                'stitch-delta-neutral',
                'stitch-delta-slight-worsening',
                'stitch-delta-strong-worsening',
                'stitch-delta-failed'
            ];

            line.classList.remove(...classes);

            const previousTitle = line.querySelector(
                'title[data-stitch-difference]'
            );

            previousTitle?.remove();

            if (mode !== 'difference') return;

            const id = String(
                line.dataset.edgeId ||
                line.dataset.failedEdgeId ||
                ''
            );

            const baselineEdge = baselineById.get(id);
            const currentEdge = currentById.get(id);

            if (!baselineEdge) return;

            const failed = !currentEdge ||
                (window.stitchTopologyFailed || [])
                    .map(String).includes(id);

            let classification = 'stitch-delta-neutral';
            let description = 'No measurable change';

            if (failed) {
                classification = 'stitch-delta-failed';
                description = 'Link failed / unavailable';
            } else {
                const before = Number(baselineEdge.utilization);
                const after = Number(currentEdge.utilization);

                if (
                    Number.isFinite(before) &&
                    Number.isFinite(after) &&
                    before > 0
                ) {
                    const delta = 100 * (after - before) / before;

                    if (delta <= -20) {
                        classification = 'stitch-delta-strong-improvement';
                    } else if (delta < -1) {
                        classification = 'stitch-delta-improvement';
                    } else if (delta <= 1) {
                        classification = 'stitch-delta-neutral';
                    } else if (delta <= 20) {
                        classification = 'stitch-delta-slight-worsening';
                    } else {
                        classification = 'stitch-delta-strong-worsening';
                    }

                    description =
                        `${delta >= 0 ? '+' : ''}${delta.toFixed(2)}%` +
                        ` utilization change` +
                        ` (${before.toFixed(3)} → ${after.toFixed(3)})`;
                } else {
                    description = 'Relative change unavailable';
                }
            }

            line.classList.add(classification);
        });

        $('stitch-topology-canvas-legend').textContent =
            mode === 'difference'
                ? 'DARK GREEN: HIGH IMPROVEMENT · GREEN: IMPROVED · GRAY: UNCHANGED · YELLOW: SLIGHTLY WORSE · RED: MUCH WORSE · DASHED: FAILED'
                : mode === 'compare'
                    ? 'LEFT: BASELINE · RIGHT: CURRENT SCENARIO'
                    : '● VERTICES · — PHYSICAL LINKS';

        const hasBaseline = baseline.length > 0;

        $('stitch-stat-delta').textContent =
            !hasBaseline
                ? 'Baseline unavailable'
                : mode === 'single'
                    ? 'Select Difference for edge deltas'
                    : 'Measured edge differences';
    }

    function update() {
        const all = links(), svg = $('graph');

        if (!svg) return;

        const vertices = new Set();

        all.forEach(l => {
            vertices.add(l.source);
            vertices.add(l.target);
        });

        $('stitch-topology-count').textContent =
            `${vertices.size} NODES · ${all.length} EDGES`;

        const over = all.filter(
            l => num(l.utilization) !== null &&
                num(l.utilization) > 1
        );

        const max = all.reduce(
            (m, l) => Math.max(m, num(l.utilization) || 0),
            0
        );

        $('stitch-stat-max').textContent =
            all.length ? max.toFixed(3) + '×' : '—';

        $('stitch-stat-over').textContent = String(over.length);

        $('stitch-stat-selected').textContent =
            String(window.stitchGetSelectedIds?.().length || 0);

        $('stitch-topology-state').textContent =
            mode === 'single'
                ? 'STATE: CURRENT TOPOLOGY'
                : mode === 'compare'
                    ? 'STATE: BASELINE VS CURRENT'
                    : 'STATE: DIFFERENTIAL Δ';

        $('stitch-topology-canvas-label').textContent =
            mode === 'single'
                ? 'VIEW: SINGLE TOPOLOGY'
                : mode === 'compare'
                    ? 'VIEW: BASELINE / CURRENT'
                    : 'VIEW: DIFFERENCE Δ (CURRENT − BASELINE)';

        const byId = new Map(
            all.map(l => [String(l.edge_id), l])
        );

        svg.querySelectorAll('line').forEach(line => {
            const id =
                line.dataset.edgeId ||
                line.dataset.failedEdgeId;

            if (!id) return;

            const l = byId.get(id);

            if (!l) return;

            const util = num(l.utilization);

            const selected =
                (window.stitchGetSelectedIds?.() || [])
                    .map(String).includes(id);

            const visible =
                filter === 'all' ||
                (
                    filter === 'congested' &&
                    util !== null &&
                    util > 1
                ) ||
                (
                    filter === 'modified' &&
                    (l.topology_failed || l.topology_unmeasured)
                ) ||
                (
                    filter === 'selected' &&
                    selected
                );

            line.classList.toggle(
                'stitch-filtered-out',
                !visible
            );

            line.style.opacity = visible ? '1' : '.08';

            line.style.stroke =
                l.topology_failed
                    ? '#ba1a1a'
                    : selected
                        ? '#005236'
                        : util !== null && util > 1
                            ? '#ba1a1a'
                            : '#a6a4a9';

            line.style.strokeWidth =
                selected
                    ? '4'
                    : l.topology_failed
                        ? '3'
                        : '2';

            line.style.strokeDasharray =
                l.topology_failed ? '7 4' : '';
        });

        paintMode(svg, all);

        const body = $('stitch-audit-rows');

        body.replaceChildren();

        all.slice(0, 1000).forEach(l => {
            const tr = document.createElement('tr');

            [
                l.edge_id,
                l.source,
                l.target,
                l.capacity,
                num(l.utilization) === null
                    ? '—'
                    : num(l.utilization).toFixed(3),
                l.topology_failed
                    ? 'FAILED'
                    : num(l.utilization) > 1
                        ? 'OVERLOADED'
                        : 'ACTIVE'
            ].forEach(v => {
                const td = document.createElement('td');

                td.textContent = String(v ?? '—');

                tr.append(td);
            });

            tr.addEventListener('click', () =>
                svg.querySelector(
                    `line[data-edge-id="${CSS.escape(String(l.edge_id))}"],line[data-failed-edge-id="${CSS.escape(String(l.edge_id))}"]`
                )?.dispatchEvent(
                    new MouseEvent('click', {bubbles: true})
                )
            );

            body.append(tr);
        });
    }

    // Preserve the existing public refresh interface.
    window.stitchRefreshTopology = update;

    document.addEventListener('DOMContentLoaded', () => {
        if (!document.body.classList.contains('topology-view')) {
            const g = $('graph');
            const wrap = document.querySelector('.graph-panel .graph-wrap');

            if (g && wrap && g.parentElement !== wrap) {
                wrap.appendChild(g);
            }
        }

        document.querySelectorAll(
            '#stitch-topology-ribbon [data-mode]'
        ).forEach(b =>
            b.addEventListener('click', () => {
                mode = b.dataset.mode;

                document.querySelectorAll(
                    '#stitch-topology-ribbon [data-mode]'
                ).forEach(x =>
                    x.classList.toggle('active', x === b)
                );

                update();
            })
        );

        $('stitch-topology-filter')?.addEventListener('change', e => {
            filter = e.target.value;
            update();
        });

        $('stitch-topology-rerun')?.addEventListener(
            'click',
            () => $('nav-analysis')?.click()
        );

        $('stitch-telemetry-toggle')?.addEventListener('click', () => {
            const el = $('stitch-telemetry-body');

            el.hidden = !el.hidden;
            el.style.display = el.hidden ? 'none' : 'grid';
        });

        $('stitch-audit-toggle')?.addEventListener('click', () => {
            const el = $('stitch-topology-audit');

            el.classList.toggle('expanded');

            $('stitch-audit-toggle').textContent =
                el.classList.contains('expanded')
                    ? 'COLLAPSE DETAILS'
                    : 'EXPAND DETAILS';
        });

        $('stitch-topology-audit-header')?.addEventListener(
            'click',
            e => {
                if (e.target.id !== 'stitch-audit-toggle') {
                    $('stitch-audit-toggle').click();
                }
            }
        );

        $('graph')?.addEventListener(
            'click',
            () => requestAnimationFrame(update)
        );

        update();
    });
})();
