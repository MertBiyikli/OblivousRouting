async function executeScenario(request, networkId, comparison = false) {
    const endpoint = comparison ? "comparisons" : "scenarios";
    const url = `/api/v1/networks/${encodeURIComponent(networkId)}/${endpoint}`;

    const response = await fetch(url, {
        method: "POST", headers: {"Content-Type": "application/json"}, body: JSON.stringify(request)
    });

    const data = await response.json();

    if (!response.ok) {
        throw new Error(typeof data.detail === "string" ? data.detail : `API returned HTTP ${response.status}`);
    }

    return data;
}