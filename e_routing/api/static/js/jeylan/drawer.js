function setJeylanMode(mode) {
    if (mode !== 'ask' && mode !== 'experiment') {
        throw new Error('Unknown Jeylan interaction mode');
    }

    const ask = mode === 'ask';

    $('jeylan-read-panel').hidden = !ask;
    $('jeylan-experiment-panel').hidden = ask;

    $('jeylan-mode-ask').classList.toggle(
        'active', ask
    );

    $('jeylan-mode-experiment').classList.toggle(
        'active', !ask
    );

    $('jeylan-mode-ask').setAttribute(
        'aria-pressed', String(ask)
    );

    $('jeylan-mode-experiment').setAttribute(
        'aria-pressed', String(!ask)
    );
}

function setJeylanOpen(open, {focus = false} = {}) {
    const isOpen = Boolean(open);

    document.body.classList.toggle('jeylan-open', isOpen);
    document.body.classList.toggle('jeylan-minimized', !isOpen);

    const panel = $('jeylan-panel');
    const launcher = $('jeylan-launcher');
    const closeButton = $('topology-jeylan-toggle');

    panel.inert = !isOpen;

    launcher.setAttribute('aria-expanded', String(isOpen));
    closeButton.setAttribute('aria-expanded', String(isOpen));

    if (focus) {
        if (isOpen) {
            closeButton.focus();
        } else {
            launcher.focus();
        }
    }
}