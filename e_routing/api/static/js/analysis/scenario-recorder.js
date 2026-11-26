/** Scenario recording and session run history. */
(function () {
    let history = [], lastEvidence = '', pending = null, previousRunText = '';
    const $ = id => document.getElementById(id);

    function record(status, message) {
        history.unshift({status, message, time: new Date().toLocaleTimeString()});
        document.dispatchEvent(new CustomEvent('jupiter:scenario-recorded', {
            detail: {
                id: document.getElementById('scenario-id')?.value || '',
                status, type: document.getElementById('analysis-type')?.value || '',
                name: message, network: document.getElementById('network')?.value || '',
                timestamp: new Date().toISOString(),
                baseline: document.getElementById('baseline-congestion')?.textContent || '',
                scenario: document.getElementById('scenario-congestion')?.textContent || '',
                delta: document.getElementById('capacity-improvement')?.textContent || '',
                evidence: document.getElementById('evidence')?.textContent || ''
            }
        }));
        history = history.slice(0, 40);
        $('stitch-history-count').textContent = history.length + ' RUN' + (history.length === 1 ? '' : 'S');
        $('stitch-history-items').replaceChildren();
        history.forEach(item => {
            const row = document.createElement('div');
            const label = document.createElement('span');
            label.textContent = item.time + ' · ' + item.message;
            const state = document.createElement('span');
            state.textContent = item.status.toUpperCase();
            state.className = item.status === 'completed' ? 'stitch-history-ok' : 'stitch-history-error';
            row.append(label, state);
            $('stitch-history-items').append(row);
        });
    }

    document.addEventListener('DOMContentLoaded', () => {
        const run = $('run'), ev = $('evidence'), status = $('status');
        run?.addEventListener('click', () => {
            pending = {
                type: $('analysis-type')?.selectedOptions?.[0]?.textContent?.trim() || 'Scenario',
                at: Date.now()
            };
            previousRunText = ev?.textContent?.trim() || '';
        });
        // Evidence is updated by the existing engine response handlers. Record only
        // substantive updates following a run, not simulated or preset log entries.
        if (ev) new MutationObserver(() => {
            if (!pending) return;
            const message = ev.textContent.trim();
            if (!message || message === previousRunText || message === lastEvidence ||
                /^(No |Waiting|Configure |Run a scenario)/i.test(message)) return;
            const currentStatus = status?.textContent?.trim() || '';
            if (/error|failed/i.test(currentStatus)) {
                record('failed', pending.type);
                lastEvidence = message;
                pending = null;
                return;
            }
            if (/running|loading|executing/i.test(currentStatus)) return;
            record('completed', pending.type);
            lastEvidence = message;
            pending = null;
        }).observe(ev, {childList: true, characterData: true, subtree: true});
        const errors = $('error');
        if (errors) new MutationObserver(() => {
            if (!pending) return;
            const message = errors.textContent.trim();
            if (message && !errors.hidden) {
                record('failed', pending.type);
                pending = null;
            }
        }).observe(errors, {
            childList: true,
            characterData: true,
            subtree: true,
            attributes: true,
            attributeFilter: ['hidden']
        });
    });
})();