function setDisplayedResultProvenance({
                                          origin, scenarioId, edgeId = null, capacityFactor = null
                                      }) {
    const element = $('displayed-result-provenance');

    if (!scenarioId) {
        clearDisplayedResultProvenance();
        return;
    }

    if (origin === 'jeylan') {
        element.textContent = `Viewing Jeylan experiment · ` + `Edge ${edgeId} ×${capacityFactor} · ` + `Workspace unchanged · ` + `Analysis ${scenarioId}`;

    } else {
        element.textContent = `Viewing Scenario Studio analysis · ` + `${scenarioId}`;
    }

    element.hidden = false;
}


function clearDisplayedResultProvenance() {
    const element = $('displayed-result-provenance');
    element.hidden = true;
    element.textContent = '';
}