// Scenario Studio: call renderScenarioAnalysis once after the scenario POST completes.
// Assumes an existing <div id="scenario-analysis"></div> in the Studio page.
const fmt = (number, digits = 2) => typeof number === 'number' && Number.isFinite(number)
  ? number.toFixed(digits) : 'N/A';
const pct = (fraction) => typeof fraction === 'number' && Number.isFinite(fraction)
  ? `${(Math.abs(fraction) < 1e-12 ? 0 : fraction * 100).toFixed(2)}%` : 'N/A';

async function renderScenarioAnalysis(networkId, scenarioId) {
  const panel = document.getElementById('scenario-analysis');
  if (!panel) throw new Error('Missing #scenario-analysis container');
  panel.replaceChildren();
  const url = `/api/v1/networks/${encodeURIComponent(networkId)}/scenarios/${encodeURIComponent(scenarioId)}/analysis`;
  const response = await fetch(url);
  if (!response.ok) {
    panel.textContent = response.status === 404
      ? 'Scenario analysis has not been generated.'
      : `Unable to load scenario analysis (HTTP ${response.status}).`;
    return;
  }
  const report = await response.json();
  const h2 = document.createElement('h2');
  h2.textContent = 'Capacity Upgrade Analysis';
  panel.append(h2);
  for (const solver of report.solver_results) {
    const section = document.createElement('section');
    const h3 = document.createElement('h3');
    h3.textContent = solver.solver_algorithms.join(', ');
    section.append(h3);
    const m = solver.performance_metrics;
    const summary = document.createElement('p');
    summary.textContent = `Congestion: ${fmt(m.baseline_congestion)} → ${fmt(m.fixed_routing_congestion)} (fixed) → ${fmt(m.reoptimized_congestion)} (reoptimized). Capacity-only improvement: ${pct(m.capacity_only_improvement)}. Additional rerouting gain: ${pct(m.reoptimization_gain)}.`;
    section.append(summary);
    const extra = document.createElement('p');
    extra.textContent = `Changed physical links: ${m.changed_physical_links ?? 'N/A'}. Total absolute load redistribution: ${fmt(m.total_absolute_load_redistribution)}.`;
    section.append(extra);
    const h4 = document.createElement('h4');
    h4.textContent = 'Reoptimized bottlenecks';
    section.append(h4);
    const ul = document.createElement('ul');
    for (const edge of solver.bottlenecks.reoptimized) {
      const li = document.createElement('li');
      li.textContent = `Edge ${edge.edge_id}: load ${fmt(edge.load)}, capacity ${fmt(edge.capacity)}, utilization ${fmt(edge.utilization)}`;
      ul.append(li);
    }
    section.append(ul);
    panel.append(section);
  }
}
