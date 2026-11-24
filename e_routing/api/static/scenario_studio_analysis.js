
const fmt = (value, digits = 2) =>
    typeof value === "number" && Number.isFinite(value)
        ? value.toFixed(digits)
        : "N/A";

function analysisElement(tag, content, parent) {
  const element = document.createElement(tag);
  element.textContent = content;
  parent.appendChild(element);
  return element;
}

async function renderScenarioAnalysis(networkId, scenarioId) {
  const panel = document.getElementById("scenario-analysis");

  if (!panel) {
    throw new Error("Missing scenario-analysis container");
  }

  panel.replaceChildren();

  const url =
      `/api/v1/networks/${encodeURIComponent(networkId)}` +
      `/scenarios/${encodeURIComponent(scenarioId)}/analysis`;

  const response = await fetch(url);

  if (!response.ok) {
    throw new Error(
        `Unable to load analysis (HTTP ${response.status})`
    );
  }

  const report = await response.json();

  const isFailure = report.scenario_kind === "link_failure";
  const title = isFailure
      ? "Link Failure & Resilience Analysis"
      : "Capacity Planning Analysis";

  analysisElement("h2", title, panel);

  if (isFailure) {
    analysisElement(
        "p",
        `Failed physical links: ${report.link_failures.join(", ")}`,
        panel
    );
  } else {
    const upgrades = report.capacity_upgrades.map(
        upgrade =>
            `Edge ${upgrade.edge_id} ×${fmt(upgrade.capacity_factor)}`
    );

    analysisElement(
        "p",
        `Capacity upgrades: ${upgrades.join(", ")}`,
        panel
    );
  }

  for (const solver of report.solver_results) {
    const section = document.createElement("section");
    panel.appendChild(section);

    analysisElement(
        "h3",
        solver.solver_algorithms.join(", "),
        section
    );

    const baseline = solver.baseline;
    const fixed = solver.fixed_routing;
    const modified = solver.reoptimized;

    analysisElement(
        "p",
        `Baseline congestion: ${fmt(baseline.max_congestion)}`,
        section
    );

    if (fixed) {
      analysisElement(
          "p",
          `Fixed-routing congestion: ${fmt(fixed.max_congestion)}`,
          section
      );
    }

    analysisElement(
        "p",
        `Reoptimized congestion: ${fmt(modified.max_congestion)}`,
        section
    );

    const before = baseline.max_congestion;
    const after = modified.max_congestion;

    const absoluteChange = after - before;
    const relativeChange = before > 0
        ? 100 * absoluteChange / before
        : null;

    analysisElement(
        "p",
        `Congestion change: ${absoluteChange >= 0 ? "+" : ""}` +
        `${fmt(absoluteChange)} ` +
        `(${relativeChange !== null && relativeChange >= 0 ? "+" : ""}` +
        `${fmt(relativeChange)}%)`,
        section
    );

    analysisElement(
        "h4",
        "Most congested links after modification",
        section
    );

    const list = document.createElement("ul");

    for (const link of solver.bottlenecks.reoptimized) {
      analysisElement(
          "li",
          `Edge ${link.edge_id}: ` +
          `load ${fmt(link.load)}, ` +
          `capacity ${fmt(link.capacity)}, ` +
          `utilization ${fmt(link.utilization)}`,
          list
      );
    }

    if (list.children.length > 4) {
      const overflow = document.createElement('details');
      const summary = document.createElement('summary');
      summary.textContent = `Show remaining ${list.children.length - 4} congested links`;
      const more = document.createElement('ul');
      while (list.children.length > 4) more.appendChild(list.children[4]);
      overflow.append(summary, more);
      section.append(list, overflow);
    } else {
      section.appendChild(list);
    }

    if (isFailure) {
      analysisElement(
          "h4",
          "Traffic redistribution",
          section
      );

      const deltas = [...solver.link_deltas]
          .sort(
              (a, b) =>
                  Math.abs(b.delta_from_baseline) -
                  Math.abs(a.delta_from_baseline)
          )
          .slice(0, 5);

      const deltaList = document.createElement("ul");

      for (const link of deltas) {
        analysisElement(
            "li",
            `Edge ${link.edge_id}: ` +
            `${link.delta_from_baseline >= 0 ? "+" : ""}` +
            fmt(link.delta_from_baseline) +
            " load units",
            deltaList
        );
      }

      section.appendChild(deltaList);
    }
  }
}


// F3-H: bounded N-1 outcome browsing. Every outcome remains accessible.
function renderN1Analysis(report, onSelect) {
  const panel = document.getElementById('scenario-analysis');
  panel.hidden = false;
  panel.replaceChildren();
  analysisElement('h2', 'N−1 Network Resilience', panel);

  const summary = report.summary;
  analysisElement('p',
      `${summary.tested_links} physical links evaluated · ` +
      `${summary.successful} successful recomputations · ` +
      `${summary.disconnected} disconnections · ` +
      `${summary.execution_failed} execution failures`, panel);
  analysisElement('p',
      `Baseline maximum congestion: ${fmt(report.baseline.max_congestion, 3)}`,
      panel);

  const baseline = report.baseline.max_congestion;
  const allOutcomes = [...report.outcomes].sort((a, b) => {
    const priority = {disconnected: 0, success: 1, execution_failed: 2};
    const difference = (priority[a.status] ?? 3) - (priority[b.status] ?? 3);
    if (difference) return difference;
    if (a.status === 'success')
      return (b.reoptimized?.max_congestion ?? -Infinity) -
          (a.reoptimized?.max_congestion ?? -Infinity);
    return a.edge_id - b.edge_id;
  });

  const tools = document.createElement('div');
  tools.className = 'result-browser-toolbar';
  const filterLabel = analysisElement('label', 'Filter failure outcomes', tools);
  const filter = document.createElement('select');
  filter.setAttribute('aria-label', 'Filter N-1 outcomes');
  for (const [value, title] of [
    ['all', 'All outcomes'],
    ['disconnected', 'Disconnected'],
    ['success', 'Successful recomputations'],
    ['execution_failed', 'Execution failures']
  ]) {
    const option = document.createElement('option');
    option.value = value;
    option.textContent = title;
    filter.appendChild(option);
  }
  filterLabel.htmlFor = 'n1-filter';
  filter.id = 'n1-filter';
  tools.appendChild(filter);
  const count = analysisElement('span', '', tools);
  count.className = 'result-browser-count';
  panel.appendChild(tools);

  const list = document.createElement('div');
  list.id = 'n1-failure-ranking';
  list.className = 'bounded-result-list';
  panel.appendChild(list);

  const pager = document.createElement('div');
  pager.className = 'result-browser-pager';
  const previous = document.createElement('button');
  const next = document.createElement('button');
  const pageInfo = document.createElement('span');
  previous.type = next.type = 'button';
  previous.className = next.className = 'secondary result-page-button';
  previous.textContent = 'Previous';
  next.textContent = 'Next';
  const expand = document.createElement('button');
  expand.type = 'button';
  expand.className = 'secondary result-expand-button';
  pager.append(previous, pageInfo, next, expand);
  panel.appendChild(pager);

  let expanded = false;
  let page = 0;
  let selectedEdge = null;
  const PAGE_SIZE = 20;
  const COLLAPSED_SIZE = 4;

  function filteredOutcomes() {
    return filter.value === 'all' ? allOutcomes :
        allOutcomes.filter(item => item.status === filter.value);
  }

  function selectOutcome(outcome) {
    selectedEdge = outcome.edge_id;
    onSelect(outcome);
    list.querySelectorAll('button[data-n1-edge-id]').forEach(button => {
      button.setAttribute('aria-pressed',
          String(Number(button.dataset.n1EdgeId) === selectedEdge));
    });
  }

  function renderPage() {
    const items = filteredOutcomes();
    const pageSize = expanded ? PAGE_SIZE : COLLAPSED_SIZE;
    const pages = Math.max(1, Math.ceil(items.length / pageSize));
    page = Math.min(page, pages - 1);
    const from = page * pageSize;
    list.replaceChildren();
    for (const outcome of items.slice(from, from + pageSize)) {
      const button = document.createElement('button');
      button.type = 'button';
      button.className = 'bottleneck-row';
      button.dataset.n1EdgeId = String(outcome.edge_id);
      button.setAttribute('aria-pressed', String(outcome.edge_id === selectedEdge));
      const label = analysisElement('span',
          `Edge ${outcome.edge_id} · ${outcome.source} → ${outcome.target}`,
          button);
      const value = document.createElement('strong');
      if (outcome.status === 'disconnected') value.textContent = 'Disconnected';
      else if (outcome.status === 'execution_failed') value.textContent = 'Analysis failed';
      else {
        const after = outcome.reoptimized.max_congestion;
        const change = baseline > 0 ? 100 * (after - baseline) / baseline : null;
        value.textContent = `${fmt(after, 2)} · ` +
            (change === null ? 'N/A' : `${change >= 0 ? '+' : ''}${fmt(change)}%`);
      }
      button.appendChild(value);
      button.addEventListener('click', () => selectOutcome(outcome));
      list.appendChild(button);
    }
    count.textContent = `${items.length} matching · showing ${items.length ? from + 1 : 0}–${Math.min(from + pageSize, items.length)}`;
    pageInfo.textContent = `Page ${page + 1} of ${pages}`;
    previous.disabled = page === 0;
    next.disabled = page >= pages - 1;
    expand.textContent = expanded ? 'Show four per page' : 'Show 20 per page';
    // If the selected result is on a different page, keep its graph selection;
    // changing the browsing page must not change the displayed analysis.
  }

  previous.addEventListener('click', () => { page--; renderPage(); });
  next.addEventListener('click', () => { page++; renderPage(); });
  expand.addEventListener('click', () => { expanded = !expanded; page = 0; renderPage(); });
  filter.addEventListener('change', () => { page = 0; renderPage(); });
  renderPage();
  if (allOutcomes.length) selectOutcome(allOutcomes[0]);
}
