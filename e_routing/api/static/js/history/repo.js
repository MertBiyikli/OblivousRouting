/**
 * History Repository
 *
 * Responsibilities:
 * - Maintain session-local scenario records.
 * - Filter, sort and render recorded scenarios.
 * - Display scenario metrics and execution evidence.
 * - Export individual records or complete session history.
 * - Integrate History with application navigation.
 *
 * Consumes:
 * - jupiter:scenario-recorded
 * - Existing History DOM elements.
 * - Existing topology SVG.
 */

(function () {
    const $ = id => document.getElementById(id);
    const rows = [];
    let selected = null, type = 'ALL', query = '', status = 'ALL', descending = true;

    const types = {
        capacity_upgrade: 'CAPACITY', link_failure: 'FAILURE', n1_resilience: 'N1', traffic_growth: 'TRAFFIC'
    };

    const names = {
        CAPACITY: 'Capacity Upgrade', FAILURE: 'Link Failure', N1: 'N-1 Resilience', TRAFFIC: 'Traffic Growth'
    };

    function go(view) {
        $('nav-' + view)?.click();

        document.body.classList.toggle('stitch-history-view', view === 'history');

        if (view === 'history') {
            document.body.classList.remove('topology-view');
            $('nav-topology')?.classList.remove('active');
            $('nav-analysis')?.classList.remove('active');
            $('nav-history')?.classList.add('active');
        } else {
            $('nav-history')?.classList.remove('active');
        }
    }

    $('nav-history').addEventListener('click', () => go('history'));

    ['topology', 'analysis'].forEach(v => $('nav-' + v)?.addEventListener('click', () => {
        document.body.classList.remove('stitch-history-view');
        $('nav-history').classList.remove('active');
    }));

    $('sh-new').addEventListener('click', () => go('analysis'));

    function download(name, obj) {
        const blob = new Blob([JSON.stringify(obj, null, 2)], {type: 'application/json'});

        const u = URL.createObjectURL(blob);
        const a = document.createElement('a');

        a.href = u;
        a.download = name;
        a.click();

        setTimeout(() => URL.revokeObjectURL(u), 1000);
    }

    function fmt(v) {
        return v && v.trim() !== '—' ? v.trim() : '—';
    }

    function render() {
        $('sh-count').textContent = rows.length + ' SCENARIOS RECORDED';

        const filtered = rows
            .filter(r => (type === 'ALL' || r.classification === type) && (status === 'ALL' || r.status === status) && [r.name, r.id, r.network]
                .join(' ')
                .toLowerCase()
                .includes(query))
            .sort((a, b) => descending ? b.time - a.time : a.time - b.time);

        $('sh-rows').replaceChildren();

        $('sh-empty').hidden = filtered.length > 0;

        $('sh-pagination').textContent = 'SHOWING ' + filtered.length + ' OF ' + rows.length + ' SCENARIOS';

        filtered.forEach((r, i) => {
            const tr = document.createElement('tr');

            if (selected === r) {
                tr.className = 'selected';
            }

            const values = [String(i + 1), r.name, names[r.classification] || r.classification, r.network || '—', new Date(r.time).toLocaleString(), r.delta || '—', r.status.toUpperCase()];

            values.forEach((v, j) => {
                const td = document.createElement('td');

                if (j === 1) {
                    td.textContent = v;

                    const sm = document.createElement('small');
                    sm.textContent = r.id || 'Session run';

                    td.append(sm);
                } else if (j === 2) {
                    const tag = document.createElement('span');

                    tag.className = 'sh-tag';
                    tag.textContent = v;

                    td.append(tag);
                } else {
                    td.textContent = v;
                }

                if (j === 6) {
                    td.className = r.status === 'completed' ? 'sh-ok' : 'sh-bad';
                }

                tr.append(td);
            });

            const op = document.createElement('td');
            const btn = document.createElement('button');

            btn.textContent = 'EXPORT';

            btn.addEventListener('click', e => {
                e.stopPropagation();
                download((r.id || 'scenario') + '.json', r);
            });

            op.append(btn);
            tr.append(op);

            tr.addEventListener('click', () => {
                selected = r;
                inspector();
                render();
            });

            $('sh-rows').append(tr);
        });
    }

    function inspector() {
        if (!selected) return;

        const r = selected;

        $('sh-inspector-type').textContent = names[r.classification] || r.classification;

        $('sh-inspector-status').textContent = r.status.toUpperCase();

        $('sh-inspector-name').textContent = r.name;

        $('sh-inspector-time').textContent = new Date(r.time).toLocaleString() + ' · ' + (r.id || 'Session run');

        $('sh-baseline').textContent = fmt(r.baseline);
        $('sh-scenario').textContent = fmt(r.scenario);
        $('sh-delta').textContent = fmt(r.delta);
        $('sh-network').textContent = r.network || '—';

        $('sh-log').textContent = r.evidence || 'No execution evidence was captured.';

        const mini = $('sh-mini-graph');

        mini.replaceChildren();

        // A topology preview is only valid for the
        // currently displayed scenario.
        if (r === rows[0] && $('graph')?.querySelectorAll('line').length) {
            const svg = $('graph').cloneNode(true);

            svg.removeAttribute('id');

            svg.querySelectorAll('.stitch-commodity-overlay,.stitch-delta-label').forEach(e => e.remove());

            mini.append(svg);
        } else {
            mini.textContent = 'Topology snapshot not stored for this run.';
        }
    }

    document.addEventListener('jupiter:scenario-recorded', e => {
        const d = e.detail || {};

        const rec = {
            ...d, classification: types[d.type] || 'OTHER', time: Date.parse(d.timestamp) || Date.now()
        };

        rows.unshift(rec);

        if (rows.length > 100) {
            rows.pop();
        }

        selected = rec;

        render();
        inspector();
    });

    $('sh-search').addEventListener('input', e => {
        query = e.target.value.toLowerCase().trim();
        render();
    });

    $('sh-types').addEventListener('click', e => {
        const b = e.target.closest('[data-type]');

        if (!b) return;

        type = b.dataset.type;

        $('sh-types').querySelectorAll('button')
            .forEach(x => x.classList.toggle('selected', x === b));

        render();
    });

    $('sh-state').addEventListener('change', e => {
        status = e.target.value;
        render();
    });

    $('sh-sort').addEventListener('click', () => {
        descending = !descending;

        $('sh-sort').textContent = descending ? 'TIMESTAMP ↓' : 'TIMESTAMP ↑';

        render();
    });

    $('sh-export-all').addEventListener('click', () => download('jupiter-session-history.json', rows));

    $('sh-export-one').addEventListener('click', () => {
        if (selected) {
            download((selected.id || 'scenario') + '.json', selected);
        }
    });

    $('sh-open-topology').addEventListener('click', () => go('topology'));

    $('sh-close').addEventListener('click', () => {
        selected = null;

        $('sh-inspector-name').textContent = 'Select a scenario';

        $('sh-mini-graph').textContent = 'Select a recorded scenario to inspect its topology.';

        render();
    });

    render();
})();
