function setGraphView(x, y, width, height) {
    const minWidth = 90, maxWidth = 1440;
    width = Math.min(maxWidth, Math.max(minWidth, width));
    height = width * GRAPH_BOUNDS.height / GRAPH_BOUNDS.width;
    graphView = {x, y, width, height};
    $('graph').setAttribute('viewBox', `${x} ${y} ${width} ${height}`);
    const baselineGraph = document.getElementById('stitch-baseline-graph');
    if (baselineGraph) baselineGraph.setAttribute('viewBox', `${x} ${y} ${width} ${height}`);
}

function resetGraphView() {
    const svg = $('graph');
    const nodes = [...svg.querySelectorAll('circle')]
        .filter(node =>
            !node.classList.contains('stitch-commodity-overlay')
        );

    if (!nodes.length) {
        setGraphView(0, 0, GRAPH_BOUNDS.width, GRAPH_BOUNDS.height);
        return;
    }

    const xs = nodes.map(node => Number(node.getAttribute('cx')));
    const ys = nodes.map(node => Number(node.getAttribute('cy')));

    const minX = Math.min(...xs);
    const maxX = Math.max(...xs);
    const minY = Math.min(...ys);
    const maxY = Math.max(...ys);

    const svgRect = svg.getBoundingClientRect();
    const inspector = $('topology-inspector');

    let reservedRight = 0;

    if (document.body.classList.contains('topology-view') &&
        inspector) {
        const inspectorRect = inspector.getBoundingClientRect();

        if (inspectorRect.width > 0 &&
            inspectorRect.right > svgRect.left &&
            inspectorRect.left < svgRect.right) {
            reservedRight = Math.max(
                0,
                svgRect.right - inspectorRect.left + 16
            );
        }
    }

    const usableWidth = Math.max(
        1,
        svgRect.width - reservedRight
    );

    const padding = 30;
    const contentWidth = Math.max(1, maxX - minX + 2 * padding);
    const contentHeight = Math.max(1, maxY - minY + 2 * padding);

    const aspect = GRAPH_BOUNDS.width / GRAPH_BOUNDS.height;

    const width = Math.max(
        contentWidth,
        contentHeight * aspect
    );

    const height = width / aspect;

    const centerX = (minX + maxX) / 2;
    const centerY = (minY + maxY) / 2;

    const shiftX = reservedRight > 0
        ? (reservedRight / Math.max(1, svgRect.width)) * width / 2
        : 0;

    setGraphView(
        centerX - width / 2 + shiftX,
        centerY - height / 2,
        width,
        height
    );
}

function zoomGraph(factor, anchorX = 360, anchorY = 170) {
    const old = graphView;
    const width = Math.min(1440, Math.max(90, old.width * factor));
    const ratio = width / old.width;
    setGraphView(
        anchorX - (anchorX - old.x) * ratio,
        anchorY - (anchorY - old.y) * ratio,
        width, old.height * ratio
    );
}

function graphCoordinates(svg, clientX, clientY) {
    const point = svg.createSVGPoint();
    point.x = clientX;
    point.y = clientY;

    const transform = svg.getScreenCTM();

    if (!transform) return null;

    return point.matrixTransform(transform.inverse());
}

function attachGraphNavigation(svg) {
    if (!svg) return;

    svg.style.touchAction = 'none';

    svg.addEventListener('wheel', event => {
        event.preventDefault();

        const anchor = graphCoordinates(
            svg,
            event.clientX,
            event.clientY
        );

        if (!anchor) return;

        zoomGraph(
            event.deltaY < 0 ? 0.86 : 1.16,
            anchor.x,
            anchor.y
        );
    }, {passive: false});

    svg.addEventListener('pointerdown', event => {
        if (event.button !== 0 || graphPan) return;

        // Keep link and vertex selection available.
        if (event.target.closest('line, circle, text')) return;

        const anchor = graphCoordinates(
            svg,
            event.clientX,
            event.clientY
        );

        if (!anchor) return;

        graphPan = {
            svg,
            pointerId: event.pointerId,
            anchor
        };

        svg.setPointerCapture(event.pointerId);
    });

    svg.addEventListener('pointermove', event => {
        if (!graphPan ||
            graphPan.svg !== svg ||
            graphPan.pointerId !== event.pointerId) {
            return;
        }

        const current = graphCoordinates(
            svg,
            event.clientX,
            event.clientY
        );

        if (!current) return;

        setGraphView(
            graphView.x + graphPan.anchor.x - current.x,
            graphView.y + graphPan.anchor.y - current.y,
            graphView.width,
            graphView.height
        );
    });

    function stopPan(event) {
        if (!graphPan ||
            graphPan.svg !== svg ||
            graphPan.pointerId !== event.pointerId) {
            return;
        }

        graphPan = null;

        if (svg.hasPointerCapture(event.pointerId)) {
            svg.releasePointerCapture(event.pointerId);
        }
    }

    svg.addEventListener('pointerup', stopPan);
    svg.addEventListener('pointercancel', stopPan);
    svg.addEventListener('lostpointercapture', stopPan);
}