(() => {
    "use strict";
    const legendMinimum =
        document.querySelector(
            "#legend-min"
        );

    const legendMiddle =
        document.querySelector(
            "#legend-mid"
        );

    const legendMaximum =
        document.querySelector(
            "#legend-max"
        );

    const fileInput =
        document.querySelector("#file-input");

    const candidateFileInput =
        document.querySelector(
            "#candidate-file-input"
        );

    const visualizationModeSelect =
        document.querySelector(
            "#visualization-mode"
        );

    const comparisonDescription =
        document.querySelector(
            "#comparison-description"
        );

    const comparisonMetrics =
        document.querySelector(
            "#comparison-metrics"
        );

    const comparisonRankings =
        document.querySelector(
            "#comparison-rankings"
        );

    const svg =
        document.querySelector("#network");

    const emptyState =
        document.querySelector("#empty-state");

    const fitButton =
        document.querySelector("#fit-button");

    const resetButton =
        document.querySelector("#reset-button");

    const toggleEdgeButton =
        document.querySelector("#toggle-edge");

    const colorModeSelect =
        document.querySelector(
            "#color-mode"
        );

    const colorScaleDescription =
        document.querySelector(
            "#color-scale-description"
        );

    const edgeCapacityInput =
        document.querySelector(
            "#edge-capacity-input"
        );

    const applyCapacityButton =
        document.querySelector(
            "#apply-capacity"
        );

    const increaseCapacityButton =
        document.querySelector(
            "#increase-capacity"
        );

    const doubleCapacityButton =
        document.querySelector(
            "#double-capacity"
        );

    const resetCapacityButton =
        document.querySelector(
            "#reset-capacity"
        );

    const capacityChangeDescription =
        document.querySelector(
            "#capacity-change-description"
        );


    let data = null;
    let selectedEdge = null;
    let candidateData = null;

    let visualizationMode =
        "utilization";

    let comparisonScale = {
        maximumAbsoluteDelta: 1
    };

    let positions =
        new Map();

    let drag = null;

    let colorMode = "relative";

    let colorScale = {
        lower: 0,
        upper: 1
    };

    let scenarioMode = "baseline";

    let view = {
        x: 0,
        y: 0,
        width: 1000,
        height: 700
    };

    function hasScenarioChanges() {
        if (!data) {
            return false;
        }

        return data.edges.some(
            edge =>
                edge.enabled === false ||
                edge.modified === true
        );
    }

    function applyEdgeCapacity(
        edge,
        newCapacity
    ) {
        const capacity =
            Number(newCapacity);

        if (
            !Number.isFinite(capacity) ||
            capacity <= 0
        ) {
            alert(
                "Capacity must be a positive number."
            );

            return;
        }

        edge.capacity =
            capacity;

        /*
         * Static scenario:
         * routing and load remain unchanged.
         */
        edge.load =
            edge.originalLoad;

        edge.utilization =
            edge.load / edge.capacity;

        edge.modified =
            Math.abs(
                edge.capacity -
                edge.originalCapacity
            ) > 1e-12;

        scenarioMode =
            hasScenarioChanges()
                ? "modified"
                : "baseline";

        updateColorScale();
        updateLegend();
        updatePanels();
        selectEdge(edge);
    }

    function resetEdgeCapacity(edge) {
        edge.capacity =
            edge.originalCapacity;

        edge.load =
            edge.originalLoad;

        edge.utilization =
            edge.originalUtilization;

        edge.modified =
            false;

        scenarioMode =
            hasScenarioChanges()
                ? "modified"
                : "baseline";

        updateColorScale();
        updateLegend();
        updatePanels();
        selectEdge(edge);
    }

    function physicalEdgeKey(edge) {
        const source =
            String(edge.source);

        const target =
            String(edge.target);

        return source < target
            ? `${source}::${target}`
            : `${target}::${source}`;
    }

    function validateComparableResults(
        baseline,
        candidate
    ) {
        if (
            baseline.nodes.length !==
            candidate.nodes.length
        ) {
            throw new Error(
                "Baseline and candidate contain different numbers of nodes."
            );
        }

        if (
            baseline.edges.length !==
            candidate.edges.length
        ) {
            throw new Error(
                "Baseline and candidate contain different numbers of physical edges."
            );
        }

        const baselineNodes =
            new Set(
                baseline.nodes.map(node =>
                    String(node.id)
                )
            );

        for (const node of candidate.nodes) {
            if (
                !baselineNodes.has(
                    String(node.id)
                )
            ) {
                throw new Error(
                    `Candidate contains node ${node.id}, which is absent from the baseline.`
                );
            }
        }

        const baselineEdges =
            new Set(
                baseline.edges.map(
                    physicalEdgeKey
                )
            );

        for (const edge of candidate.edges) {
            const key =
                physicalEdgeKey(edge);

            if (!baselineEdges.has(key)) {
                throw new Error(
                    `Candidate contains edge ${edge.source} ↔ ${edge.target}, which is absent from the baseline.`
                );
            }
        }
    }

    function buildCandidateEdgeIndex() {
        const index =
            new Map();

        if (!candidateData) {
            return index;
        }

        candidateData.edges.forEach(edge => {
            index.set(
                physicalEdgeKey(edge),
                edge
            );
        });

        return index;
    }

    function calculateComparison() {
        if (!data || !candidateData) {
            return;
        }

        const candidateEdges =
            buildCandidateEdgeIndex();

        let maximumAbsoluteDelta = 0;

        data.edges.forEach(edge => {
            const candidateEdge =
                candidateEdges.get(
                    physicalEdgeKey(edge)
                );

            if (!candidateEdge) {
                edge.comparison = null;
                return;
            }

            const baselineLoad =
                Number(edge.originalLoad ?? edge.load);

            const candidateLoad =
                Number(candidateEdge.load);

            const baselineUtilization =
                Number(
                    edge.originalUtilization ??
                    edge.utilization
                );

            const candidateUtilization =
                Number(candidateEdge.utilization);

            const loadDelta =
                candidateLoad -
                baselineLoad;

            const utilizationDelta =
                candidateUtilization -
                baselineUtilization;

            const relativeLoadDelta =
                Math.abs(baselineLoad) > 1e-12
                    ? loadDelta /
                    Math.abs(baselineLoad)
                    : (
                        Math.abs(candidateLoad) > 1e-12
                            ? 1
                            : 0
                    );

            edge.comparison = {
                baselineLoad,
                candidateLoad,
                loadDelta,
                relativeLoadDelta,

                baselineCapacity:
                    Number(
                        edge.originalCapacity ??
                        edge.capacity
                    ),

                candidateCapacity:
                    Number(candidateEdge.capacity),

                baselineUtilization,
                candidateUtilization,
                utilizationDelta,

                candidateEdgeId:
                candidateEdge.id
            };

            maximumAbsoluteDelta =
                Math.max(
                    maximumAbsoluteDelta,
                    Math.abs(utilizationDelta)
                );
        });

        comparisonScale.maximumAbsoluteDelta =
            maximumAbsoluteDelta > 1e-12
                ? maximumAbsoluteDelta
                : 1;
    }

    function validateResult(result) {
        if (
            result === null ||
            typeof result !== "object"
        ) {
            throw new Error(
                "The selected file does not contain a JSON object."
            );
        }

        if (!Array.isArray(result.nodes)) {
            throw new Error(
                "The visualization result does not contain a nodes array."
            );
        }

        if (!Array.isArray(result.edges)) {
            throw new Error(
                "The visualization result does not contain an edges array."
            );
        }

        if (
            result.summary === null ||
            typeof result.summary !== "object"
        ) {
            throw new Error(
                "The visualization result does not contain summary metrics."
            );
        }
    }

    function buildLayout() {
        positions = new Map();

        if (!data || data.nodes.length === 0) {
            return;
        }

        const centerX = 500;
        const centerY = 350;

        const radius =
            Math.min(
                280,
                100 + data.nodes.length * 7
            );

        data.nodes.forEach((node, index) => {
            const angle =
                (
                    2 *
                    Math.PI *
                    index
                ) /
                data.nodes.length -
                Math.PI / 2;

            positions.set(
                String(node.id),
                {
                    x:
                        centerX +
                        Math.cos(angle) * radius,

                    y:
                        centerY +
                        Math.sin(angle) * radius
                }
            );
        });

        /*
         * A lightweight force-relaxation pass.
         *
         * This is not intended to replace Cytoscape or Graphviz, but it makes
         * small and medium topologies more readable without dependencies.
         */
        const iterations = 180;

        for (
            let iteration = 0;
            iteration < iterations;
            ++iteration
        ) {
            const forces =
                new Map();

            data.nodes.forEach(node => {
                forces.set(
                    String(node.id),
                    { x: 0, y: 0 }
                );
            });

            /*
             * Node repulsion.
             */
            for (
                let firstIndex = 0;
                firstIndex < data.nodes.length;
                ++firstIndex
            ) {
                for (
                    let secondIndex = firstIndex + 1;
                    secondIndex < data.nodes.length;
                    ++secondIndex
                ) {
                    const firstNode =
                        data.nodes[firstIndex];

                    const secondNode =
                        data.nodes[secondIndex];

                    const first =
                        positions.get(
                            String(firstNode.id)
                        );

                    const second =
                        positions.get(
                            String(secondNode.id)
                        );

                    let deltaX =
                        first.x - second.x;

                    let deltaY =
                        first.y - second.y;

                    let distanceSquared =
                        deltaX * deltaX +
                        deltaY * deltaY;

                    if (distanceSquared < 1) {
                        deltaX =
                            Math.random() - 0.5;

                        deltaY =
                            Math.random() - 0.5;

                        distanceSquared = 1;
                    }

                    const distance =
                        Math.sqrt(distanceSquared);

                    const repulsion =
                        9000 / distanceSquared;

                    const forceX =
                        repulsion *
                        deltaX /
                        distance;

                    const forceY =
                        repulsion *
                        deltaY /
                        distance;

                    const firstForce =
                        forces.get(
                            String(firstNode.id)
                        );

                    const secondForce =
                        forces.get(
                            String(secondNode.id)
                        );

                    firstForce.x += forceX;
                    firstForce.y += forceY;

                    secondForce.x -= forceX;
                    secondForce.y -= forceY;
                }
            }

            /*
             * Edge attraction.
             */
            data.edges.forEach(edge => {
                const source =
                    positions.get(
                        String(edge.source)
                    );

                const target =
                    positions.get(
                        String(edge.target)
                    );

                if (!source || !target) {
                    return;
                }

                const deltaX =
                    target.x - source.x;

                const deltaY =
                    target.y - source.y;

                const distance =
                    Math.max(
                        1,
                        Math.sqrt(
                            deltaX * deltaX +
                            deltaY * deltaY
                        )
                    );

                const desiredLength = 130;

                const attraction =
                    (
                        distance -
                        desiredLength
                    ) * 0.012;

                const forceX =
                    attraction *
                    deltaX /
                    distance;

                const forceY =
                    attraction *
                    deltaY /
                    distance;

                const sourceForce =
                    forces.get(
                        String(edge.source)
                    );

                const targetForce =
                    forces.get(
                        String(edge.target)
                    );

                if (sourceForce) {
                    sourceForce.x += forceX;
                    sourceForce.y += forceY;
                }

                if (targetForce) {
                    targetForce.x -= forceX;
                    targetForce.y -= forceY;
                }
            });

            /*
             * Pull everything gently toward the canvas centre.
             */
            data.nodes.forEach(node => {
                const key =
                    String(node.id);

                const position =
                    positions.get(key);

                const force =
                    forces.get(key);

                force.x +=
                    (centerX - position.x) *
                    0.002;

                force.y +=
                    (centerY - position.y) *
                    0.002;

                const cooling =
                    1 -
                    iteration / iterations;

                position.x +=
                    Math.max(
                        -12,
                        Math.min(
                            12,
                            force.x * cooling
                        )
                    );

                position.y +=
                    Math.max(
                        -12,
                        Math.min(
                            12,
                            force.y * cooling
                        )
                    );
            });
        }
    }

    function findFailureAnalysis(edge) {
        if (
            !data ||
            !Array.isArray(data.linkFailures)
        ) {
            return null;
        }

        return data.linkFailures.find(
            failure =>
                Number(failure.failedEdgeId) ===
                Number(edge.id) ||
                Number(failure.antiEdgeId) ===
                Number(edge.id)
        ) ?? null;
    }
    function percentile(sortedValues, percentileValue) {
        if (sortedValues.length === 0) {
            return 0;
        }

        const position =
            (sortedValues.length - 1) *
            percentileValue;

        const lowerIndex =
            Math.floor(position);

        const upperIndex =
            Math.ceil(position);

        if (lowerIndex === upperIndex) {
            return sortedValues[lowerIndex];
        }

        const fraction =
            position - lowerIndex;

        return (
            sortedValues[lowerIndex] *
            (1 - fraction)
        ) + (
            sortedValues[upperIndex] *
            fraction
        );
    }

    function updateLegend() {
        if (
            visualizationMode === "difference" &&
            candidateData
        ) {
            legendMinimum.textContent =
                "Improved";

            legendMiddle.textContent =
                "Unchanged";

            legendMaximum.textContent =
                "Worse";

            return;
        }
        if (colorMode === "capacity") {
            legendMinimum.textContent =
                "0%";

            legendMiddle.textContent =
                "70%";

            legendMaximum.textContent =
                "100%+";

            return;
        }

        const middle =
            (
                colorScale.lower +
                colorScale.upper
            ) / 2;

        legendMinimum.textContent =
            formatPercent(
                colorScale.lower
            );

        legendMiddle.textContent =
            formatPercent(middle);

        legendMaximum.textContent =
            formatPercent(
                colorScale.upper
            );
    }

    function updateColorScale() {
        if (!data || data.edges.length === 0) {
            colorScale = {
                lower: 0,
                upper: 1
            };

            return;
        }

        const utilizations =
            data.edges
                /*
                 * A failed edge should not influence the scenario color scale.
                 */
                .filter(edge =>
                    edge.enabled !== false
                )
                .map(edge =>
                    Number(edge.utilization)
                )
                .filter(value =>
                    Number.isFinite(value) &&
                    value >= 0
                )
                .sort(
                    (first, second) =>
                        first - second
                );

        if (utilizations.length === 0) {
            colorScale = {
                lower: 0,
                upper: 1
            };

            return;
        }

        const lower =
            percentile(
                utilizations,
                0.05
            );

        const upper =
            percentile(
                utilizations,
                0.95
            );

        colorScale = {
            lower,
            upper:
                upper > lower
                    ? upper
                    : lower + 1
        };
    }
    function clamp(value, minimum, maximum) {
        return Math.max(
            minimum,
            Math.min(maximum, value)
        );
    }

    function edgeColor(utilization) {
        const value =
            Number(utilization);

        if (!Number.isFinite(value)) {
            return "hsl(0 0% 60%)";
        }

        /*
         * Absolute capacity-oriented interpretation.
         */
        if (colorMode === "capacity") {
            if (value < 0.7) {
                return "hsl(142 70% 38%)";
            }

            if (value < 0.85) {
                return "hsl(48 90% 48%)";
            }

            if (value < 1.0) {
                return "hsl(25 95% 52%)";
            }

            return "hsl(0 72% 48%)";
        }

        /*
         * Relative percentile-normalized interpretation.
         *
         * Apply log1p because congestion values can differ by several orders
         * of magnitude.
         */
        const transformedValue =
            Math.log1p(
                Math.max(0, value)
            );

        const transformedLower =
            Math.log1p(
                Math.max(
                    0,
                    colorScale.lower
                )
            );

        const transformedUpper =
            Math.log1p(
                Math.max(
                    0,
                    colorScale.upper
                )
            );

        const denominator =
            transformedUpper -
            transformedLower;

        const normalized =
            denominator > 0
                ? clamp(
                    (
                        transformedValue -
                        transformedLower
                    ) / denominator,
                    0,
                    1
                )
                : 0.5;

        /*
         * Hue:
         *
         * 120 = green
         * 60  = yellow
         * 30  = orange
         * 0   = red
         */
        const hue =
            120 * (1 - normalized);

        return `hsl(${hue} 78% 45%)`;
    }

    function comparisonEdgeColor(edge) {
        if (!edge.comparison) {
            return "hsl(0 0% 60%)";
        }

        const delta =
            Number(
                edge.comparison.utilizationDelta
            );

        if (!Number.isFinite(delta)) {
            return "hsl(0 0% 60%)";
        }

        const normalizedMagnitude =
            clamp(
                Math.abs(delta) /
                comparisonScale.maximumAbsoluteDelta,
                0,
                1
            );

        /*
         * Candidate improved the edge.
         */
        if (delta < -1e-12) {
            const lightness =
                72 -
                normalizedMagnitude * 34;

            return `hsl(142 68% ${lightness}%)`;
        }

        /*
         * Candidate degraded the edge.
         */
        if (delta > 1e-12) {
            const lightness =
                72 -
                normalizedMagnitude * 34;

            return `hsl(0 72% ${lightness}%)`;
        }

        /*
         * No meaningful change.
         */
        return "hsl(215 10% 62%)";
    }

    function formatPercent(value) {
        if (!Number.isFinite(Number(value))) {
            return "—";
        }

        return `${(
            Number(value) * 100
        ).toFixed(1)}%`;
    }

    function formatSignedPercent(value) {
        const numericValue =
            Number(value);

        if (!Number.isFinite(numericValue)) {
            return "—";
        }

        const sign =
            numericValue > 0
                ? "+"
                : "";

        return `${sign}${(
            numericValue * 100
        ).toFixed(1)}%`;
    }

    function formatSignedNumber(value) {
        const numericValue =
            Number(value);

        if (!Number.isFinite(numericValue)) {
            return "—";
        }

        const sign =
            numericValue > 0
                ? "+"
                : "";

        return `${sign}${formatNumber(
            numericValue
        )}`;
    }

    function formatNumber(value) {
        if (!Number.isFinite(Number(value))) {
            return "—";
        }

        return Number(value).toLocaleString(
            undefined,
            {
                maximumFractionDigits: 3
            }
        );
    }

    function render() {
        svg.replaceChildren();

        svg.setAttribute(
            "viewBox",
            [
                view.x,
                view.y,
                view.width,
                view.height
            ].join(" ")
        );

        if (!data) {
            return;
        }

        const edgesGroup =
            document.createElementNS(
                "http://www.w3.org/2000/svg",
                "g"
            );

        const nodesGroup =
            document.createElementNS(
                "http://www.w3.org/2000/svg",
                "g"
            );

        data.edges.forEach(edge => {
            const source =
                positions.get(
                    String(edge.source)
                );

            const target =
                positions.get(
                    String(edge.target)
                );

            if (!source || !target) {
                return;
            }

            /*
             * Wide transparent line gives the edge a larger click target.
             */
            const hit =
                document.createElementNS(
                    "http://www.w3.org/2000/svg",
                    "line"
                );

            setLine(hit, source, target);

            hit.classList.add("edge-hit");

            hit.addEventListener(
                "click",
                () => selectEdge(edge)
            );

            edgesGroup.appendChild(hit);

            const line =
                document.createElementNS(
                    "http://www.w3.org/2000/svg",
                    "line"
                );

            setLine(line, source, target);

            line.classList.add("edge-line");

            if (
                selectedEdge &&
                selectedEdge.id === edge.id
            ) {
                line.classList.add("selected");
            }

            line.style.stroke =
                edge.enabled === false
                    ? "#94a3b8"
                    : (
                        visualizationMode === "difference" &&
                        candidateData
                            ? comparisonEdgeColor(edge)
                            : edgeColor(edge.utilization)
                    );

            line.style.strokeWidth =
                String(
                    2.5 +
                    Math.min(
                        3,
                        Math.sqrt(
                            Math.max(
                                0,
                                Number(edge.capacity)
                            )
                        ) * 0.15
                    )
                );

            line.style.opacity =
                edge.enabled === false
                    ? "0.22"
                    : "0.82";

            line.style.strokeDasharray =
                edge.enabled === false
                    ? "9 7"
                    : "";

            line.addEventListener(
                "click",
                () => selectEdge(edge)
            );

            edgesGroup.appendChild(line);

            if (edge.modified) {
                line.classList.add("modified");
            }
            if (edge.enabled === false) {
                line.style.strokeDasharray =
                    "9 7";
            } else if (edge.modified) {
                line.style.strokeDasharray =
                    "5 3";
            } else {
                line.style.strokeDasharray = "";
            }
            });

        data.nodes.forEach(node => {
            const position =
                positions.get(
                    String(node.id)
                );

            if (!position) {
                return;
            }

            const circle =
                document.createElementNS(
                    "http://www.w3.org/2000/svg",
                    "circle"
                );

            circle.setAttribute(
                "cx",
                position.x
            );

            circle.setAttribute(
                "cy",
                position.y
            );

            circle.setAttribute(
                "r",
                "16"
            );

            circle.classList.add("node");

            nodesGroup.appendChild(circle);

            const label =
                document.createElementNS(
                    "http://www.w3.org/2000/svg",
                    "text"
                );

            label.setAttribute(
                "x",
                position.x
            );

            label.setAttribute(
                "y",
                position.y
            );

            label.textContent =
                node.label ?? node.id;

            label.classList.add(
                "node-label"
            );

            nodesGroup.appendChild(label);
        });

        svg.append(
            edgesGroup,
            nodesGroup
        );
    }

    function setLine(
        line,
        source,
        target
    ) {
        line.setAttribute(
            "x1",
            source.x
        );

        line.setAttribute(
            "y1",
            source.y
        );

        line.setAttribute(
            "x2",
            target.x
        );

        line.setAttribute(
            "y2",
            target.y
        );
    }

    function updateComparisonPanels() {
        if (!candidateData) {
            comparisonMetrics.hidden =
                true;

            comparisonRankings.hidden =
                true;

            comparisonDescription.textContent =
                "Load a candidate result to enable comparison.";

            return;
        }

        comparisonMetrics.hidden =
            false;

        comparisonRankings.hidden =
            false;

        const summary =
            calculateComparisonSummary();

        document.querySelector(
            "#comparison-baseline-maximum"
        ).textContent =
            formatPercent(
                summary.baselineMaximum
            );

        document.querySelector(
            "#comparison-candidate-maximum"
        ).textContent =
            formatPercent(
                summary.candidateMaximum
            );

        document.querySelector(
            "#comparison-maximum-change"
        ).textContent =
            formatSignedPercent(
                summary.maximumDelta
            );

        document.querySelector(
            "#comparison-average-change"
        ).textContent =
            formatSignedPercent(
                summary.averageDelta
            );

        const baselineSolver =
            data.solver || "Baseline";

        const candidateSolver =
            candidateData.solver || "Candidate";

        comparisonDescription.textContent =
            `${baselineSolver} compared with ${candidateSolver}. Green links improved; red links became worse.`;

        const comparableEdges =
            data.edges.filter(edge =>
                edge.comparison !== null &&
                edge.comparison !== undefined
            );

        const improved =
            [...comparableEdges]
                .filter(edge =>
                    edge.comparison.utilizationDelta <
                    -1e-12
                )
                .sort(
                    (first, second) =>
                        first.comparison.utilizationDelta -
                        second.comparison.utilizationDelta
                )
                .slice(0, 8);

        const degraded =
            [...comparableEdges]
                .filter(edge =>
                    edge.comparison.utilizationDelta >
                    1e-12
                )
                .sort(
                    (first, second) =>
                        second.comparison.utilizationDelta -
                        first.comparison.utilizationDelta
                )
                .slice(0, 8);

        updateComparisonRanking(
            "#improved-edge-ranking",
            improved,
            "No links improved."
        );

        updateComparisonRanking(
            "#degraded-edge-ranking",
            degraded,
            "No links degraded."
        );
    }

    function updatePanels() {
        const scenarioDetails =
            document.querySelector(
                "#scenario-details"
            );

        scenarioDetails.innerHTML = `
            <div>
                <dt>Graph</dt>
                <dd>${escapeHtml(
            data.graphName || "Unnamed"
        )}</dd>
            </div>

            <div>
                <dt>Solver</dt>
                <dd>${escapeHtml(
            data.solver || "Unknown"
        )}</dd>
            </div>

            <div>
                <dt>Demand</dt>
                <dd>${escapeHtml(
            data.demandModel || "Unknown"
        )}</dd>
            </div>
        `;

        const metrics =
            document.querySelectorAll(
                "#metrics strong"
            );

        const summary =
            calculateScenarioSummary();

        metrics[0].textContent =
            formatPercent(
                summary.maximumCongestion
            );

        metrics[1].textContent =
            formatPercent(
                summary.averageUtilization
            );

        metrics[2].textContent =
            formatNumber(
                summary.overloadedEdges
            );

        metrics[3].textContent =
            formatNumber(
                summary.totalRoutedDemand
            );

        const baselineMaximum =
            document.querySelector(
                "#baseline-maximum"
            );

        const baselineAverage =
            document.querySelector(
                "#baseline-average"
            );

        const baselineOverloaded =
            document.querySelector(
                "#baseline-overloaded"
            );

        if (baselineMaximum) {
            baselineMaximum.textContent =
                `Baseline ${formatPercent(
                    data.summary.maximumCongestion
                )}`;
        }

        if (baselineAverage) {
            baselineAverage.textContent =
                `Baseline ${formatPercent(
                    data.summary.averageUtilization
                )}`;
        }

        if (baselineOverloaded) {
            baselineOverloaded.textContent =
                `Baseline ${formatNumber(
                    data.summary.overloadedEdges
                )}`;
        }

        const ranking =
            [...data.edges]
                .sort(
                    (first, second) =>
                        second.utilization -
                        first.utilization
                )
                .slice(0, 8);

        const rankingElement =
            document.querySelector(
                "#edge-ranking"
            );

        rankingElement.innerHTML =
            ranking.map(edge => `
                <div
                    class="ranking-row utilization-ranking-row"
                    data-edge-id="${edge.id}"
                >
                    <span>
                        ${edge.source}
                        ↔
                        ${edge.target}
                    </span>

                    <span>
                        ${formatPercent(
                edge.utilization
            )}
                    </span>
                </div>
            `).join("") || "No edges.";

        document
            .querySelectorAll(".utilization-ranking-row")
            .forEach(row => {
                row.addEventListener(
                    "click",
                    () => {
                        const edge =
                            data.edges.find(
                                item =>
                                    String(item.id) ===
                                    row.dataset.edgeId
                            );

                        if (edge) {
                            selectEdge(edge);
                        }
                    }
                );
            });

        const failureRankingElement =
            document.querySelector(
                "#failure-ranking"
            );

        if (
            Array.isArray(data.linkFailures) &&
            data.linkFailures.length > 0
        ) {
            const failures =
                [...data.linkFailures]
                    .sort(
                        (first, second) =>
                            second.lostTraffic -
                            first.lostTraffic
                    )
                    .slice(0, 8);

            failureRankingElement.innerHTML =
                failures.map(failure => `
            <div
                class="ranking-row failure-ranking-row"
                data-edge-id="${failure.failedEdgeId}"
            >
                <span>
                    ${failure.source}
                    ↔
                    ${failure.target}
                </span>

                <span>
                    ${formatNumber(
                    failure.lostTraffic
                )}
                </span>
            </div>
        `).join("");

            document
                .querySelectorAll(
                    ".failure-ranking-row"
                )
                .forEach(row => {
                    row.addEventListener(
                        "click",
                        () => {
                            const edge =
                                data.edges.find(
                                    item =>
                                        Number(item.id) ===
                                        Number(
                                            row.dataset.edgeId
                                        )
                                );

                            if (edge) {
                                selectEdge(edge);
                            }
                        }
                    );
                });
        } else {
            failureRankingElement.textContent =
                "No failure analysis available.";
        }

        const capacityRankingElement =
            document.querySelector(
                "#capacity-ranking"
            );

        if (capacityRankingElement) {
            const targetUtilization =
                0.8;

            const recommendations =
                data.edges
                    .filter(edge =>
                        edge.enabled !== false
                    )
                    .map(edge => {
                        const currentCapacity =
                            Number(edge.capacity);

                        const recommendedCapacity =
                            requiredCapacity(
                                edge,
                                targetUtilization
                            );

                        const capacityIncrease =
                            Math.max(
                                0,
                                recommendedCapacity -
                                currentCapacity
                            );

                        const relativeIncrease =
                            currentCapacity > 0
                                ? capacityIncrease /
                                currentCapacity
                                : 0;

                        return {
                            edge,
                            recommendedCapacity,
                            capacityIncrease,
                            relativeIncrease
                        };
                    })
                    .filter(recommendation =>
                        recommendation.capacityIncrease >
                        1e-12
                    )
                    .sort(
                        (first, second) =>
                            second.relativeIncrease -
                            first.relativeIncrease
                    )
                    .slice(0, 8);

            if (recommendations.length === 0) {
                capacityRankingElement.textContent =
                    "No upgrades required for an 80% utilization target.";
            } else {
                capacityRankingElement.innerHTML =
                    recommendations
                        .map(recommendation => `
                    <div
                        class="ranking-row capacity-ranking-row"
                        data-edge-id="${recommendation.edge.id}"
                    >
                        <span>
                            ${recommendation.edge.source}
                            ↔
                            ${recommendation.edge.target}
                        </span>

                        <span>
                            ${formatNumber(
                            recommendation.recommendedCapacity
                        )}
                            <small>
                                +${formatPercent(
                            recommendation.relativeIncrease
                        )}
                            </small>
                        </span>
                    </div>
                `)
                        .join("");

                document
                    .querySelectorAll(
                        ".capacity-ranking-row"
                    )
                    .forEach(row => {
                        row.addEventListener(
                            "click",
                            () => {
                                selectEdgeById(
                                    row.dataset.edgeId
                                );
                            }
                        );
                    });
            }
        }

        updateComparisonPanels();
    }

    function selectEdge(edge) {
        selectedEdge = edge;

        document.querySelector(
            "#edge-empty"
        ).hidden = true;

        document.querySelector(
            "#edge-details"
        ).hidden = false;

        document.querySelector(
            "#edge-name"
        ).textContent =
            `${edge.source} ↔ ${edge.target}`;

        document.querySelector(
            "#edge-capacity"
        ).textContent =
            formatNumber(edge.capacity);

        document.querySelector(
            "#edge-load"
        ).textContent =
            formatNumber(edge.load);

        document.querySelector(
            "#edge-utilization"
        ).textContent =
            formatPercent(edge.utilization);

        document.querySelector(
            "#edge-id"
        ).textContent =
            edge.id;

        const status =
            document.querySelector(
                "#edge-status"
            );

        status.textContent =
            edge.enabled === false
                ? "Hidden"
                : "Active";

        status.classList.toggle(
            "hidden",
            edge.enabled === false
        );

        toggleEdgeButton.textContent =
            edge.enabled === false
                ? "Restore edge"
                : "Simulate edge failure";

        toggleEdgeButton.classList.toggle(
            "restore",
            edge.enabled === false
        );

        edgeCapacityInput.value =
            Number(edge.capacity);

        const capacityDifference =
            Number(edge.capacity) -
            Number(edge.originalCapacity);

        const relativeDifference =
            edge.originalCapacity > 0
                ? capacityDifference /
                edge.originalCapacity
                : 0;

        if (edge.modified) {
            const sign =
                relativeDifference >= 0
                    ? "+"
                    : "";

            capacityChangeDescription.textContent =
                `Baseline ${formatNumber(
                    edge.originalCapacity
                )}; change ${sign}${formatPercent(
                    relativeDifference
                )}.`;

            resetCapacityButton.disabled =
                false;
        } else {
            capacityChangeDescription.textContent =
                `Baseline capacity: ${formatNumber(
                    edge.originalCapacity
                )}.`;

            resetCapacityButton.disabled =
                true;
        }

        const failure =
            findFailureAnalysis(edge);

        const failureDetails =
            document.querySelector(
                "#failure-details"
            );

        if (failure) {
            failureDetails.hidden = false;

            document.querySelector(
                "#failure-commodities"
            ).textContent =
                formatNumber(
                    failure.affectedCommodities
                );

            document.querySelector(
                "#failure-demand"
            ).textContent =
                formatNumber(
                    failure.affectedDemandVolume
                );

            document.querySelector(
                "#failure-lost"
            ).textContent =
                formatNumber(
                    failure.lostTraffic
                );

            document.querySelector(
                "#failure-fraction"
            ).textContent =
                formatPercent(
                    failure.affectedDemandFraction
                );
        } else {
            failureDetails.hidden = true;
        }

        const comparisonDetails =
            document.querySelector(
                "#edge-comparison-details"
            );

        if (
            candidateData &&
            edge.comparison
        ) {
            comparisonDetails.hidden =
                false;

            document.querySelector(
                "#edge-baseline-load"
            ).textContent =
                formatNumber(
                    edge.comparison.baselineLoad
                );

            document.querySelector(
                "#edge-candidate-load"
            ).textContent =
                formatNumber(
                    edge.comparison.candidateLoad
                );

            document.querySelector(
                "#edge-load-change"
            ).textContent =
                formatSignedNumber(
                    edge.comparison.loadDelta
                );

            document.querySelector(
                "#edge-baseline-utilization"
            ).textContent =
                formatPercent(
                    edge.comparison.baselineUtilization
                );

            document.querySelector(
                "#edge-candidate-utilization"
            ).textContent =
                formatPercent(
                    edge.comparison.candidateUtilization
                );

            document.querySelector(
                "#edge-utilization-change"
            ).textContent =
                formatSignedPercent(
                    edge.comparison.utilizationDelta
                );
        } else {
            comparisonDetails.hidden =
                true;
        }
        render();
    }

    function calculateScenarioSummary() {
        if (!data) {
            return null;
        }

        const activeEdges =
            data.edges.filter(
                edge =>
                    edge.enabled !== false
            );

        let maximumUtilization = 0;
        let utilizationSum = 0;
        let overloadedEdges = 0;

        activeEdges.forEach(edge => {
            const utilization =
                Number(edge.utilization);

            if (!Number.isFinite(utilization)) {
                return;
            }

            maximumUtilization =
                Math.max(
                    maximumUtilization,
                    utilization
                );

            utilizationSum +=
                utilization;

            if (utilization > 1.0) {
                ++overloadedEdges;
            }
        });

        return {
            maximumCongestion:
            maximumUtilization,

            averageUtilization:
                activeEdges.length > 0
                    ? utilizationSum /
                    activeEdges.length
                    : 0,

            overloadedEdges,

            totalRoutedDemand:
            data.summary.totalRoutedDemand
        };
    }


    function calculateComparisonSummary() {
        if (!data || !candidateData) {
            return null;
        }

        const baselineMaximum =
            Number(
                data.summary.maximumCongestion
            );

        const candidateMaximum =
            Number(
                candidateData.summary.maximumCongestion
            );

        const baselineAverage =
            Number(
                data.summary.averageUtilization
            );

        const candidateAverage =
            Number(
                candidateData.summary.averageUtilization
            );

        return {
            baselineMaximum,
            candidateMaximum,

            maximumDelta:
                candidateMaximum -
                baselineMaximum,

            maximumRelativeChange:
                Math.abs(baselineMaximum) > 1e-12
                    ? (
                    candidateMaximum -
                    baselineMaximum
                ) / baselineMaximum
                    : 0,

            baselineAverage,
            candidateAverage,

            averageDelta:
                candidateAverage -
                baselineAverage
        };
    }

    function updateComparisonRanking(
        selector,
        edges,
        emptyMessage
    ) {
        const element =
            document.querySelector(selector);

        if (!element) {
            return;
        }

        if (edges.length === 0) {
            element.textContent =
                emptyMessage;

            return;
        }

        element.innerHTML =
            edges.map(edge => `
            <div
                class="ranking-row comparison-ranking-row"
                data-edge-id="${edge.id}"
            >
                <span>
                    ${escapeHtml(edge.source)}
                    ↔
                    ${escapeHtml(edge.target)}
                </span>

                <span>
                    ${formatSignedPercent(
                edge.comparison.utilizationDelta
            )}
                </span>
            </div>
        `).join("");

        element
            .querySelectorAll(
                ".comparison-ranking-row"
            )
            .forEach(row => {
                row.addEventListener(
                    "click",
                    () => {
                        const edge =
                            data.edges.find(
                                candidate =>
                                    String(candidate.id) ===
                                    row.dataset.edgeId
                            );

                        if (edge) {
                            selectEdge(edge);
                        }
                    }
                );
            });
    }
    function requiredCapacity(
        edge,
        targetUtilization = 0.8
    ) {
        const load =
            Number(edge.load);

        if (
            !Number.isFinite(load) ||
            load < 0 ||
            !Number.isFinite(targetUtilization) ||
            targetUtilization <= 0
        ) {
            return 0;
        }

        return load / targetUtilization;
    }

    function selectEdgeById(edgeId) {
        if (!data) {
            return;
        }

        const edge =
            data.edges.find(
                candidate =>
                    Number(candidate.id) ===
                    Number(edgeId)
            );

        if (edge) {
            selectEdge(edge);
        }
    }

    function escapeHtml(value) {
        const element =
            document.createElement("div");

        element.textContent =
            String(value);

        return element.innerHTML;
    }

    fileInput.addEventListener(
        "change",
        async event => {
            const file =
                event.target.files?.[0];

            if (!file) {
                return;
            }

            try {
                const text =
                    await file.text();

                const parsed =
                    JSON.parse(text);

                validateResult(parsed);

                data = parsed;
                candidateData =
                    null;

                visualizationMode =
                    "utilization";

                visualizationModeSelect.value =
                    "utilization";

                candidateFileInput.disabled =
                    false;
                scenarioMode =
                    "baseline";

                data.edges.forEach(edge => {
                    if (edge.enabled === undefined) {
                        edge.enabled = true;
                    }
                    edge.comparison =
                        null;

                    edge.originalCapacity =
                        Number(edge.capacity);

                    edge.originalLoad =
                        Number(edge.load);

                    edge.originalUtilization =
                        Number(edge.utilization);

                    edge.modified =
                        false;
                });

                updateColorScale();
                updateLegend();
                selectedEdge = null;

                buildLayout();

                view = {
                    x: 0,
                    y: 0,
                    width: 1000,
                    height: 700
                };

                emptyState.hidden = true;

                document.querySelector(
                    "#edge-empty"
                ).hidden = false;

                document.querySelector(
                    "#edge-details"
                ).hidden = true;

                updatePanels();
                render();
            } catch (error) {
                alert(
                    error instanceof Error
                        ? error.message
                        : String(error)
                );
            }
        }
    );

    toggleEdgeButton.addEventListener(
        "click",
        () => {
            if (!selectedEdge) {
                return;
            }

            selectedEdge.enabled =
                selectedEdge.enabled === false;

            scenarioMode =
                hasScenarioChanges()
                    ? "modified"
                    : "baseline";

            updateColorScale();
            updateLegend();
            updatePanels();
            selectEdge(selectedEdge);
        }
    );

    colorModeSelect.addEventListener(
        "change",
        event => {
            colorMode =
                event.target.value;

            updateLegend();
            render();
        }
    );

    applyCapacityButton.addEventListener(
        "click",
        () => {
            if (!selectedEdge) {
                return;
            }

            applyEdgeCapacity(
                selectedEdge,
                edgeCapacityInput.value
            );
        }
    );

    edgeCapacityInput.addEventListener(
        "keydown",
        event => {
            if (event.key !== "Enter") {
                return;
            }

            if (!selectedEdge) {
                return;
            }

            applyEdgeCapacity(
                selectedEdge,
                edgeCapacityInput.value
            );
        }
    );

    increaseCapacityButton.addEventListener(
        "click",
        () => {
            if (!selectedEdge) {
                return;
            }

            applyEdgeCapacity(
                selectedEdge,
                Number(selectedEdge.capacity) *
                1.10
            );
        }
    );

    doubleCapacityButton.addEventListener(
        "click",
        () => {
            if (!selectedEdge) {
                return;
            }

            applyEdgeCapacity(
                selectedEdge,
                Number(selectedEdge.capacity) *
                2.0
            );
        }
    );

    resetCapacityButton.addEventListener(
        "click",
        () => {
            if (!selectedEdge) {
                return;
            }

            resetEdgeCapacity(
                selectedEdge
            );
        }
    );

    resetButton.addEventListener(
        "click",
        () => {
            if (!data) {
                return;
            }

            data.edges.forEach(edge => {
                edge.enabled =
                    true;

                edge.capacity =
                    edge.originalCapacity;

                edge.load =
                    edge.originalLoad;

                edge.utilization =
                    edge.originalUtilization;

                edge.modified =
                    false;
            });

            scenarioMode =
                "baseline";

            updateColorScale();
            updateLegend();
            updatePanels();

            if (selectedEdge) {
                selectEdge(selectedEdge);
            } else {
                render();
            }
        }
    );

    candidateFileInput.addEventListener(
        "change",
        async event => {
            const file =
                event.target.files?.[0];

            if (!file) {
                return;
            }

            if (!data) {
                alert(
                    "Load a baseline result first."
                );

                candidateFileInput.value =
                    "";

                return;
            }

            try {
                const text =
                    await file.text();

                const parsed =
                    JSON.parse(text);

                validateResult(parsed);

                validateComparableResults(
                    data,
                    parsed
                );

                candidateData =
                    parsed;

                calculateComparison();

                visualizationMode =
                    "difference";

                visualizationModeSelect.value =
                    "difference";

                updateComparisonPanels();

                if (selectedEdge) {
                    selectEdge(selectedEdge);
                } else {
                    render();
                }
            } catch (error) {
                candidateData =
                    null;

                candidateFileInput.value =
                    "";

                data.edges.forEach(edge => {
                    edge.comparison =
                        null;
                });

                updateComparisonPanels();
                render();

                alert(
                    error instanceof Error
                        ? error.message
                        : String(error)
                );
            }
        }
    );

    visualizationModeSelect.addEventListener(
        "change",
        event => {
            const requestedMode =
                event.target.value;

            if (
                requestedMode === "difference" &&
                !candidateData
            ) {
                visualizationMode =
                    "utilization";

                visualizationModeSelect.value =
                    "utilization";

                alert(
                    "Load a candidate routing result before enabling difference mode."
                );

                return;
            }

            visualizationMode =
                requestedMode;

            if (
                visualizationMode === "difference"
            ) {
                comparisonDescription.textContent =
                    "Green links improved, grey links are unchanged, and red links became worse.";
            } else if (candidateData) {
                comparisonDescription.textContent =
                    "Candidate loaded. Select Candidate difference to display routing changes.";
            }

            render();
        }
    );

    fitButton.addEventListener(
        "click",
        () => {
            view = {
                x: 0,
                y: 0,
                width: 1000,
                height: 700
            };

            render();
        }
    );

    svg.addEventListener(
        "wheel",
        event => {
            if (!data) {
                return;
            }

            event.preventDefault();

            const scale =
                event.deltaY > 0
                    ? 1.12
                    : 0.89;

            const nextWidth =
                Math.min(
                    4000,
                    Math.max(
                        250,
                        view.width * scale
                    )
                );

            const nextHeight =
                nextWidth * 0.7;

            view.x +=
                (
                    view.width -
                    nextWidth
                ) / 2;

            view.y +=
                (
                    view.height -
                    nextHeight
                ) / 2;

            view.width =
                nextWidth;

            view.height =
                nextHeight;

            render();
        },
        {
            passive: false
        }
    );

    svg.addEventListener(
        "pointerdown",
        event => {
            drag = {
                x: event.clientX,
                y: event.clientY,
                viewX: view.x,
                viewY: view.y
            };

            svg.setPointerCapture(
                event.pointerId
            );
        }
    );

    svg.addEventListener(
        "pointermove",
        event => {
            if (!drag) {
                return;
            }

            view.x =
                drag.viewX -
                (
                    event.clientX -
                    drag.x
                ) *
                view.width /
                svg.clientWidth;

            view.y =
                drag.viewY -
                (
                    event.clientY -
                    drag.y
                ) *
                view.height /
                svg.clientHeight;

            render();
        }
    );

    svg.addEventListener(
        "pointerup",
        () => {
            drag = null;
        }
    );

    svg.addEventListener(
        "pointercancel",
        () => {
            drag = null;
        }
    );
})();