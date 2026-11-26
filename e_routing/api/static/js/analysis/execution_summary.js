/**
 * Execution Summary
 *
 * Presentation-only component consuming existing scenario evidence.
 *
 * Responsibilities:
 * - Display scenario execution state.
 * - Display baseline and modified maximum congestion.
 * - Display relative congestion change.
 *
 * Consumes:
 * - jupiter:scenario-recorded
 * - Existing Analysis DOM elements.
 */
(function () {
    const $ = id => document.getElementById(id);

    const text = (id, value) => {
        const el = $(id);
        if (el) el.textContent = value;
    };

    const measured = value => value !== null && value !== undefined && value !== '' && Number.isFinite(Number(value));

    const readable = value => value && value !== '—' ? value : 'N/A';

    function refreshSummary(state, detail = {}) {
        text('stitch-execution-state', state.toUpperCase());
        text('stitch-execution-type', detail.type || $('analysis-type')?.selectedOptions?.[0]?.textContent?.trim() || 'Scenario');

        if (state === 'running') {
            ['stitch-execution-baseline', 'stitch-execution-current', 'stitch-execution-delta']
                .forEach(id => text(id, '—'));

            text('stitch-execution-note', 'Awaiting engine results.');

            return;
        }

        if (state === 'failed') {
            text('stitch-execution-note', 'Execution failed. See the error message and detailed evidence in Analysis.');
            return;
        }

        const before = detail.baseline || $('baseline-congestion')?.textContent?.trim();
        const after = detail.scenario || $('scenario-congestion')?.textContent?.trim();

        text('stitch-execution-baseline', readable(before));
        text('stitch-execution-current', readable(after));

        const b = Number(before);
        const a = Number(after);

        const comparable = measured(before) && measured(after) && b > 0;

        text('stitch-execution-delta', comparable ? `${((a - b) / b * 100) >= 0 ? '+' : ''}${((a - b) / b * 100).toFixed(2)}%` : 'N/A');

        text('stitch-execution-note', comparable ? 'Change in maximum congestion relative to baseline (negative means improvement).' : 'No comparable pair of congestion measurements; see detailed scenario evidence (especially for N−1).');
    }

    document.addEventListener('DOMContentLoaded', () => {
        $('run')?.addEventListener('click', () => refreshSummary('running'));

        document.addEventListener('jupiter:scenario-recorded', event => {
            const detail = event.detail || {};

            refreshSummary(detail.status === 'completed' ? 'completed' : 'failed', detail);
        });
    });
})();
