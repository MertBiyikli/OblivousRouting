(function () {
    function syncStitchPanel() {
        const type = document.getElementById('analysis-type')?.value;
        const label = document.getElementById('stitch-current-type');
        document.body.dataset.stitchScenario = type || '';
        if (label) label.textContent = ({
            capacity_upgrade: 'CAPACITY CONFIG',
            link_failure: 'LINK FAILURE',
            n1_resilience: 'N−1 RESILIENCE',
            traffic_growth: 'TRAFFIC GROWTH'
        })[type] || 'SCENARIO CONFIG';
        const hint = document.getElementById('stitch-n1-hint');
        if (hint) hint.hidden = type !== 'n1_resilience';
    }

    document.addEventListener('DOMContentLoaded', () => {
        const select = document.getElementById('analysis-type');
        select?.addEventListener('change', syncStitchPanel);
        syncStitchPanel();
    });
})();