async function askAnalysisJeylan(context, question, selectedEdgeId, selectedEdgeIds, failedEdgeId) {
    const {networkId, analysisId, analysisType} = context;

    const url = `/api/v1/networks/${encodeURIComponent(networkId)}` + `/analyses/${encodeURIComponent(analysisType)}` + `/${encodeURIComponent(analysisId)}/ask`;

    const response = await fetch(url, {
        method: 'POST', headers: {'Content-Type': 'application/json'}, body: JSON.stringify({
            question,
            selected_edge_id: selectedEdgeId,
            selected_edge_ids: selectedEdgeIds,
            selected_failure_edge_id: failedEdgeId
        })
    });

    const result = await response.json();

    if (!response.ok) {
        throw new Error(typeof result.detail === 'string' ? result.detail : `Jeylan request failed (HTTP ${response.status})`);
    }

    if (result.analysis_id !== analysisId) {
        throw new Error('Jeylan analysis evidence mismatch');
    }

    return result;
}

async function askJeylan(question) {
    if (!activeScenarioId) {
        throw new Error("Execute a scenario before asking Jeylan.");
    }

    const requestedId = activeScenarioId;
    const response = await fetch(`/api/v1/scenarios/${encodeURIComponent(requestedId)}/ask`, {
        method: "POST", headers: {"Content-Type": "application/json"}, body: JSON.stringify({question}),
    });

    const data = await response.json();

    if (!response.ok) {
        throw new Error(typeof data.detail === "string" ? data.detail : `Jeylan request failed (${response.status})`);
    }

    return data;
}