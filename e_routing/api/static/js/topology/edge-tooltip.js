/**
 * Edge Tooltip
 *
 * Read-only presentation component for topology edge measurements.
 *
 * Responsibilities:
 * - Display capacity, load and utilization.
 * - Display baseline and current measurements.
 * - Display failed or unmeasured edge status.
 * - Position and hide the floating tooltip.
 *
 * Consumes existing topology accessors.
 */
(function () {
    const $ = id => document.getElementById(id);

    const measured = value => value !== null && value !== undefined && value !== '' && Number.isFinite(Number(value));

    const fmt = value => measured(value) ? Number(value).toFixed(3) : 'N/A';

    function tooltipText(edge, baseline, current, isFailed) {
        const lines = [`Endpoints: ${edge.source} → ${edge.target}`];

        if (isFailed) {
            lines.push('Status: FAILED / removed in scenario');
        } else if (edge.topology_unmeasured) {
            lines.push('Status: no current measurement');
        }

        const shown = isFailed ? null : edge;

        lines.push(`Capacity: ${shown && measured(shown.capacity) ? shown.capacity : 'N/A'}`);
        lines.push(`Load: ${shown && measured(shown.load) ? shown.load : 'N/A'}`);
        lines.push(`Utilization: ${shown ? fmt(shown.utilization) : 'N/A'}`);

        if (baseline) {
            lines.push(`Baseline utilization: ${fmt(baseline.utilization)}`);
        }

        if (current && !isFailed) {
            lines.push(`Current utilization: ${fmt(current.utilization)}`);
        }

        if (baseline && current && !isFailed && measured(baseline.utilization) && measured(current.utilization) && Number(baseline.utilization) > 0) {
            const delta = 100 * (Number(current.utilization) - Number(baseline.utilization)) / Number(baseline.utilization);

            lines.push(`Relative utilization Δ: ${delta >= 0 ? '+' : ''}${delta.toFixed(2)}%`);
        } else if (baseline && (!current || isFailed)) {
            lines.push('Relative utilization Δ: unavailable');
        }

        return lines.join('\n');
    }

    document.addEventListener('DOMContentLoaded', () => {
        const tooltip = document.createElement('div');

        tooltip.id = 'stitch-edge-hover';
        tooltip.setAttribute('role', 'tooltip');

        document.body.append(tooltip);

        const hide = () => {
            tooltip.style.display = 'none';
        };

        let lastLine = null;

        function move(event) {
            const line = event.target.closest?.('line[data-edge-id],line[data-failed-edge-id]');

            if (!line) {
                lastLine = null;
                hide();
                return;
            }

            const id = line.dataset.edgeId || line.dataset.failedEdgeId;

            const baseline = (window.stitchTopologyBaseline || [])
                .find(e => String(e.edge_id) === id);

            const current = (window.stitchTopologyCurrent || [])
                .find(e => String(e.edge_id) === id);

            const displayed = (window.stitchGetTopologyLinks?.() || [])
                .find(e => String(e.edge_id) === id);

            const failed = (window.stitchTopologyFailed || [])
                .some(e => String(e) === id);

            const isBaseline = event.currentTarget.id === 'stitch-baseline-graph';

            const edge = isBaseline ? baseline : (displayed || current || baseline);

            if (!edge) {
                hide();
                return;
            }

            if (line !== lastLine) {
                tooltip.replaceChildren();

                const title = document.createElement('strong');

                title.textContent = `${isBaseline ? 'BASELINE · ' : 'CURRENT · '}EDGE ${id}`;

                const info = document.createElement('span');

                info.textContent = tooltipText(edge, baseline, current, !isBaseline && failed);

                tooltip.append(title, info);

                lastLine = line;
            }

            tooltip.style.display = 'block';

            const x = Math.min(event.clientX + 16, window.innerWidth - tooltip.offsetWidth - 8);

            const y = Math.min(event.clientY + 16, window.innerHeight - tooltip.offsetHeight - 8);

            tooltip.style.left = Math.max(8, x) + 'px';
            tooltip.style.top = Math.max(8, y) + 'px';
        }

        ['graph', 'stitch-baseline-graph'].forEach(id => {
            const svg = $(id);

            if (!svg) return;

            svg.addEventListener('pointermove', move);

            svg.addEventListener('pointerleave', () => {
                lastLine = null;
                hide();
            });

            svg.addEventListener('pointerdown', hide);
        });
    });
})();
