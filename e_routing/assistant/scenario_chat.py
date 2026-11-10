
"""Evidence-grounded interpretation of a completed scenario."""

from __future__ import annotations

import json
from typing import Any

from e_routing.assistant.scenario_evidence import summarize_scenario

SCENARIO_INSTRUCTIONS = """
You are Jeylan, an evidence-grounded network planning assistant.

Answer questions only about the supplied completed scenario.
Use the scenario evidence tool before answering.
Do not invent metrics or infer unrecorded network conditions.

Congestion and utilization are load/capacity ratios unless explicitly
identified otherwise. A ratio of 1.5 corresponds to 150% utilization.

This scenario uses a persisted routing scheme and hypothetical inputs.
It is not live telemetry, a new routing optimization, or production failover.

Never claim that a scenario proves operational safety.
Cite evidence as [E1].
Respond in the user's language.
"""

ANALYSIS_INSTRUCTIONS = """
You are Jeylan, an evidence-grounded network analysis assistant.

You receive verified results from a completed E-Routing analysis.

Use the supplied evidence tool before answering.

Explain findings using only recorded metrics and outcomes.
Never invent congestion values, routing decisions, edge identities,
traffic demands, or hypothetical experiment results.

Distinguish fixed routing from reoptimized routing.
Distinguish disconnected outcomes from successful rerouting.

Congestion is a load/capacity ratio.
A ratio of 1.5 corresponds to 150% utilization.

The evidence describes modeled network behavior, not live telemetry.
Do not claim that a result proves operational safety.

If the evidence is insufficient, explain what is missing.
Cite evidence as [E1].
Respond in the user's language.

Prefer concise, structured explanations with the key result first.

Round congestion values to two decimal places for presentation.
Present relative changes as percentages, using two decimal places.
Preserve the precision of the underlying evidence.

When interpreting N-1 analysis:
- Distinguish disconnected failures from successful rerouting.
- Rank by congestion increase when asked for deterioration.
- Mention ties when multiple failures have equivalent outcomes.
- State whether changes represent increases or decreases.
- Do not interpret better post-failure routing as proof
  that removing links improves optimal network performance.

When recommending experiments:
- Describe the proposed configuration and expected measurement.
- Do not claim that an experiment was executed unless a
  validated execution tool actually returned its result.
- Distinguish observed results from hypotheses.
"""

SCENARIO_TOOL = {
    "type": "function",
    "name": "get_scenario_summary",
    "description": "Retrieve deterministic evidence for the selected completed scenario.",
    "strict": True,
    "parameters": {
        "type": "object",
        "properties": {},
        "required": [],
        "additionalProperties": False,
    },
}

def answer_evidence_question(
        evidence: dict[str, Any],
        question: str,
        client: Any,
        model: str = "gpt-5.4-mini",
        *,
        instructions: str = SCENARIO_INSTRUCTIONS,
        evidence_source: str = "persisted_routing_scenario",
        evidence_id: str | None = None,
) -> dict[str, Any]:
    if not isinstance(question, str) or not 1 <= len(question.strip()) <= 2000:
        raise ValueError("Question must contain 1–2000 characters")

    if evidence_id is None:
        evidence_id = evidence["scenario_id"]

    evidence_text = json.dumps(evidence, allow_nan=False)

    if len(evidence_text) > 50_000:
        raise ValueError("Scenario evidence exceeds the supported size")

    conversation: list[Any] = [
        {"role": "user", "content": question}
    ]
    evidence_used = False

    for _ in range(4):
        response = client.responses.create(
            model=model,
            instructions=instructions,
            input=conversation,
            tools=[SCENARIO_TOOL],
            tool_choice="required" if not evidence_used else "none",
            store=False,
        )

        if response.status != "completed":
            raise RuntimeError(
                f"Model response status: {response.status}"
            )

        conversation.extend(
            item.model_dump(exclude_none=True)
            for item in response.output
        )

        calls = [
            item for item in response.output
            if item.type == "function_call"
        ]

        if not calls:
            if not evidence_used or not response.output_text.strip():
                raise RuntimeError(
                    "Jeylan returned no evidence-backed answer"
                )

            result = {
                "answer": response.output_text,
                "evidence": [{
                    "id": "E1",
                    "source": evidence_source,
                    "scenario_id": evidence_id,
                }],
            }

            if evidence_source == "persisted_routing_scenario":
                result["scenario_id"] = evidence_id
            else:
                result["analysis_id"] = evidence_id

            return result

        for call in calls:
            if call.name != "get_scenario_summary":
                raise ValueError("Unauthorized scenario tool")

            if json.loads(call.arguments) != {}:
                raise ValueError("Unexpected scenario tool arguments")

            conversation.append({
                "type": "function_call_output",
                "call_id": call.call_id,
                "output": json.dumps({
                    "evidence_id": "E1",
                    "source": evidence_source,
                    "scenario_id": evidence_id,
                    "result": evidence,
                }, allow_nan=False),
            })

            evidence_used = True

    raise RuntimeError("Jeylan exceeded the tool-call limit")


def answer_scenario_question(
        result: dict[str, Any],
        question: str,
        client: Any,
        model: str = "gpt-5.4-mini",
) -> dict[str, Any]:
    evidence = summarize_scenario(result)

    return answer_evidence_question(
        evidence,
        question,
        client,
        model,
    )