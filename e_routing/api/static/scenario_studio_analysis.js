
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

    section.appendChild(list);

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


function renderN1Analysis(report, onSelect) {
  const panel = document.getElementById("scenario-analysis");
  panel.hidden = false;
  panel.replaceChildren();

  analysisElement("h2", "N−1 Network Resilience", panel);

  const summary = report.summary;

  analysisElement(
      "p",
      `${summary.tested_links} physical links evaluated · ` +
      `${summary.successful} successful recomputations · ` +
      `${summary.disconnected} disconnections · ` +
      `${summary.execution_failed} execution failures`,
      panel
  );

  analysisElement(
      "p",
      `Baseline maximum congestion: ${fmt(report.baseline.max_congestion, 3)}`,
      panel
  );

  const baseline = report.baseline.max_congestion;

  // Disconnections first, then the largest congestion increases.
  // Execution failures are not interpreted as network failures.
  const priority = {
    disconnected: 0,
    success: 1,
    execution_failed: 2
  };

  const outcomes = [...report.outcomes].sort((a, b) => {
    const statusDifference =
        (priority[a.status] ?? 3) - (priority[b.status] ?? 3);

    if (statusDifference) return statusDifference;

    if (a.status === "success") {
      return b.reoptimized.max_congestion -
          a.reoptimized.max_congestion;
    }

    return a.edge_id - b.edge_id;
  });

  const list = document.createElement("div");
  list.id = "n1-failure-ranking";
  list.style.display = "grid";
  list.style.gap = "8px";

  for (const outcome of outcomes) {
    const button = document.createElement("button");
    button.type = "button";
    button.className = "bottleneck-row";
    button.dataset.n1EdgeId = String(outcome.edge_id);

    const label = document.createElement("span");
    label.textContent =
        `Edge ${outcome.edge_id} · ${outcome.source} → ${outcome.target}`;

    const value = document.createElement("strong");

    if (outcome.status === "disconnected") {
      value.textContent = "Disconnected";
    } else if (outcome.status === "execution_failed") {
      value.textContent = "Analysis failed";
    } else {
      const after = outcome.reoptimized.max_congestion;
      const change = baseline > 0
          ? 100 * (after - baseline) / baseline
          : null;

      value.textContent =
          `${fmt(after, 2)} · ` +
          (change === null
              ? "N/A"
              : `${change >= 0 ? "+" : ""}${fmt(change)}%`);
    }

    button.append(label, value);
    button.addEventListener("click", () => {
      list.querySelectorAll("button").forEach(item => {
        item.setAttribute(
            "aria-pressed",
            String(item === button)
        );
      });

      onSelect(outcome);
    });

    list.appendChild(button);
  }

  panel.appendChild(list);

  // Select the highest-risk outcome initially.
  if (outcomes.length) {
    list.firstElementChild.click();
  }
}
