
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


def answer_scenario_question(
        result: dict[str, Any],
        question: str,
        client: Any,
        model: str = "gpt-5.4-mini",
) -> dict[str, Any]:
    if not isinstance(question, str) or not 1 <= len(question.strip()) <= 2000:
        raise ValueError("Question must contain 1–2000 characters")

    evidence = summarize_scenario(result)
    evidence_text = json.dumps(evidence, allow_nan=False)

    if len(evidence_text) > 50_000:
        raise ValueError("Scenario evidence exceeds the supported size")

    conversation: list[Any] = [{"role": "user", "content": question}]
    evidence_used = False

    for _ in range(4):
        response = client.responses.create(
            model=model,
            instructions=SCENARIO_INSTRUCTIONS,
            input=conversation,
            tools=[SCENARIO_TOOL],
            tool_choice="required" if not evidence_used else "none",
            store=False,
        )

        if response.status != "completed":
            raise RuntimeError(f"Model response status: {response.status}")

        conversation.extend(
            item.model_dump(exclude_none=True) for item in response.output
        )

        calls = [item for item in response.output if item.type == "function_call"]

        if not calls:
            if not evidence_used or not response.output_text.strip():
                raise RuntimeError("Jeylan returned no evidence-backed answer")

            return {
                "scenario_id": evidence["scenario_id"],
                "answer": response.output_text,
                "evidence": [
                    {
                        "id": "E1",
                        "source": "persisted_routing_scenario",
                        "scenario_id": evidence["scenario_id"],
                    }
                ],
            }

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
                    "source": "persisted_routing_scenario",
                    "scenario_id": evidence["scenario_id"],
                    "result": evidence,
                }, allow_nan=False),
            })
            evidence_used = True

    raise RuntimeError("Jeylan exceeded the tool-call limit")
