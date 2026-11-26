function svgEl(name, attrs) {
    const el = document.createElementNS('http://www.w3.org/2000/svg', name);
    for (const [k, v] of Object.entries(attrs)) el.setAttribute(k, String(v));
    return el
}

function renderBottleneckRanking(allLinks) {
    const container = $('bottleneck-list');
    container.replaceChildren();
    const selection = selectTopologyLinks(allLinks);
    if (!selection.large) return;
    const ranked = selection.ranked;
    const controls = document.createElement('div');
    controls.className = 'result-browser-pager';
    const toggle = document.createElement('button');
    toggle.type = 'button';
    toggle.className = 'secondary result-expand-button';
    const previous = document.createElement('button');
    previous.type = 'button';
    previous.textContent = 'Previous';
    previous.className = 'secondary result-page-button';
    const next = document.createElement('button');
    next.type = 'button';
    next.textContent = 'Next';
    next.className = 'secondary result-page-button';
    const pageInfo = document.createElement('span');
    controls.append(previous, pageInfo, next, toggle);
    const rows = document.createElement('div');
    rows.className = 'bounded-result-list';
    container.append(rows, controls);
    let expanded = false, page = 0;

    function render() {
        const size = expanded ? RANKING_PAGE_SIZE : INITIAL_RANKING_SIZE;
        const pages = Math.max(1, Math.ceil(ranked.length / size));
        page = Math.min(page, pages - 1);
        rows.replaceChildren();
        for (const link of ranked.slice(page * size, (page + 1) * size)) {
            const button = document.createElement('button');
            button.type = 'button';
            button.className = 'bottleneck-row';
            const label = document.createElement('span');
            label.textContent = `Edge ${link.edge_id} · ${link.source} → ${link.target}`;
            const value = document.createElement('strong');
            value.textContent = Number.isFinite(Number(link.utilization))
                ? Number(link.utilization).toFixed(3) : 'N/A';
            button.append(label, value);
            button.addEventListener('click', () => {
                const matching = [...$('graph').querySelectorAll('line')].find(
                    line => line.dataset.edgeId === String(link.edge_id));
                if (matching) matching.dispatchEvent(new Event('click'));
                else {
                    selectedLinkId = link.edge_id;
                    $('edge').value = link.edge_id;
                    $('failed-edge').value = link.edge_id;
                    $('selected').textContent = `Selected edge ${link.edge_id} (${link.source} → ${link.target})`;
                    $('ask-selected-link').disabled = !activeAnalysisContext;
                }
            });
            rows.appendChild(button);
        }
        previous.disabled = page === 0;
        next.disabled = page === pages - 1;
        pageInfo.textContent = `Page ${page + 1}/${pages} · ${ranked.length} links`;
        toggle.textContent = expanded ? 'Show four per page' : 'Show 20 per page';
    }

    previous.addEventListener('click', () => {
        page--;
        render();
    });
    next.addEventListener('click', () => {
        page++;
        render();
    });
    toggle.addEventListener('click', () => {
        expanded = !expanded;
        page = 0;
        render();
    });
    render();
}