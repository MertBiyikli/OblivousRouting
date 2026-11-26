function appendInlineMarkdown(container, value) {
    const tokens = String(value).split(/(\*\*[^*]+\*\*|`[^`]+`)/g);
    for (const token of tokens) {
        if (token.startsWith('**') && token.endsWith('**')) {
            const strong = document.createElement('strong');
            strong.textContent = token.slice(2, -2);
            container.appendChild(strong);
        } else if (token.startsWith('`') && token.endsWith('`')) {
            const code = document.createElement('code');
            code.textContent = token.slice(1, -1);
            container.appendChild(code);
        } else {
            container.appendChild(document.createTextNode(token));
        }
    }
}

function showCompositeJeylanDigest(root, result) {
    const fullBody = document.createElement('div');
    fullBody.className = 'jeylan-full-body';

    while (root.firstChild) {
        fullBody.appendChild(root.firstChild);
    }

    const digest = document.createElement('section');
    digest.className = 'jeylan-digest';

    digest.appendChild(compositeElement('strong',
        result.status === 'completed'
            ? 'Investigation completed'
            : 'Investigation partially completed'));

    digest.appendChild(compositeElement('p',
        `${result.completed_steps.length} of ${result.plan.steps.length} analyses completed. ` +
        'Simulated results; active workspace unchanged.'));

    if (result.comparison_status === 'completed') {
        for (const comparison of result.comparisons) {
            if (comparison.comparison_type !== 'capacity_n1') continue;

            const change =
                comparison.relative_baseline_congestion_change;

            const comparisonLine = compositeElement('p',
                `N−1 congestion: ${compositeNumber(comparison.before_baseline_max_congestion)} → ` +
                `${compositeNumber(comparison.after_baseline_max_congestion)} ` +
                `(${compositePercent(change)}).`);

            if (Number.isFinite(change)) {
                comparisonLine.className =
                    change < 0 ? 'favorable' :
                        change > 0 ? 'unfavorable' : '';
            }

            digest.appendChild(comparisonLine);

            const counts = comparison.counts;

            digest.appendChild(compositeElement('p',
                `${counts.improved} failure scenarios improved · ` +
                `${counts.regressed} regressed · ` +
                `${counts.became_connected} became connected · ` +
                `${counts.became_disconnected} became disconnected.`));
        }
    } else {
        digest.appendChild(compositeElement('p',
            'No verified before/after N−1 comparison available.'));
    }

    const details = document.createElement('details');
    details.className = 'jeylan-full-report';

    details.appendChild(compositeElement('summary',
        'Read full Jeylan explanation and evidence'));

    details.appendChild(fullBody);

    root.append(digest, details);
}

function renderJeylanAnswer(answer, compositeResult = null) {
    // Only round long decimal literals in the assistant's prose.
    // Persisted numerical measurements and IDs are not modified.
    const displayAnswer = String(answer).replace(
        /\b\d+\.\d{7,}\b/g,
        value => Number(value).toFixed(2)
    );
    const root = $('interpretation');
    root.replaceChildren();
    let currentList = null;
    let currentListType = null;
    for (const line of displayAnswer.split(/\r?\n/)) {
        const trimmed = line.trim();
        if (!trimmed) {
            currentList = null;
            currentListType = null;
            continue;
        }
        const heading = trimmed.match(/^#{1,4}\s+(.+)$/);
        const bullet = trimmed.match(/^[-*]\s+(.+)$/);
        const numbered = trimmed.match(/^\d+[.)]\s+(.+)$/);
        if (heading) {
            currentList = null;
            currentListType = null;
            const element = document.createElement('h3');
            appendInlineMarkdown(element, heading[1]);
            root.appendChild(element);
        } else if (bullet || numbered) {
            const kind = bullet ? 'ul' : 'ol';
            if (currentListType !== kind) {
                currentList = document.createElement(kind);
                root.appendChild(currentList);
                currentListType = kind;
            }
            const item = document.createElement('li');
            appendInlineMarkdown(item, (bullet || numbered)[1]);
            currentList.appendChild(item);
        } else {
            currentList = null;
            currentListType = null;
            const paragraph = document.createElement('p');
            appendInlineMarkdown(paragraph, trimmed);
            root.appendChild(paragraph);
        }
    }

    if (compositeResult) {
        showCompositeJeylanDigest(root, compositeResult);
    }
}