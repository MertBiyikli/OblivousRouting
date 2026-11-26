// Third-row labels match the measurements in each report, not
// the previously selected analysis type.
function setAnalysisMetricLabels(type) {
    const configurations = {
        capacity_upgrade: {
            labels: ['Baseline routing', 'Capacity-only', 'Reoptimized routing'],
            descriptions: [
                'Original capacities',
                'Fixed routing · Upgraded capacities',
                'Routing recomputed after capacity upgrade'
            ]
        },
        link_failure: {
            labels: [
                'Failed physical link',
                'Critical surviving link',
                'Largest load change'
            ],
            descriptions: [
                'Out-of-service edge',
                'Worst post-failure utilization',
                'Largest surviving-link load delta'
            ]
        },
        n1_resilience: {
            labels: [
                'Failed physical link',
                'Post-failure congestion',
                'Network impact'
            ],
            descriptions: [
                'Removed physical link',
                'N/A for disconnected or failed outcomes',
                'Disconnection or congestion change'
            ]
        },
        traffic_growth: {
            labels: ['Baseline demand', 'Modified demand', 'Demand change'],
            descriptions: [
                'Original demand volume',
                'Effective demand volume',
                'Relative change in total demand'
            ]
        }
    };

    const config = configurations[type];
    if (!config) throw new Error('Unsupported analysis type');

    for (let i = 1; i <= 3; i++) {
        $('third-label-' + i).textContent = config.labels[i - 1];
        $('third-help-' + i).textContent = config.descriptions[i - 1];
    }
}

function clearAnalysisSpecificMetrics() {
    for (const id of [
        'reference-congestion',
        'upgraded-congestion',
        'capacity-improvement'
    ]) {
        $(id).textContent = '—';
    }
    $('capacity-improvement').style.color = '';
}