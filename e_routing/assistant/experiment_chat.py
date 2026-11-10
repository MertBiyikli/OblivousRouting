
"""One bounded, evidence-grounded Jeylan capacity experiment."""

from __future__ import annotations

import json
import logging
import re
from typing import Any

from e_routing.api.workspace_models import WorkspaceSnapshot
from e_routing.assistant.bottleneck_investigation import (
    run_bottleneck_investigation,
)
from e_routing.assistant.experiment_executor import (
    run_capacity_experiment,
)
from e_routing.assistant.experiment_models import (
    CapacityExperimentPlan,
)


logger = logging.getLogger(__name__)


CAPACITY_TOOL = {
    "type": "function",
    "name": "run_capacity_experiment",
    "description": (
        "Specify one explicitly requested, simulation-only capacity "
        "upgrade. A 50 percent increase means factor 1.5. "
        "Only use this tool if the user identified a physical edge ID."
    ),
    "strict": True,
    "parameters": {
        "type": "object",
        "properties": {
            "edge_id": {
                "type": "integer",
            },
            "capacity_factor": {
                "type": "number",
            },
            "demand_model": {
                "type": "string",
                "enum": [
                    "uniform",
                    "gravity",
                    "bimodal",
                    "gaussian",
                ],
            },
        },
        "required": [
            "edge_id",
            "capacity_factor",
            "demand_model",
        ],
        "additionalProperties": False,
    },
}


BOTTLENECK_TOOL = {
    "type": "function",
    "name": "investigate_worst_bottleneck",
    "description": (
        "Find the physical link with the highest baseline "
        "utilization from a validated completed network analysis, "
        "then simulate one capacity upgrade. The backend "
        "selects the physical edge deterministically. "
        "Use when the user requests the worst bottleneck "
        "without specifying an edge ID."
    ),
    "strict": True,
    "parameters": {
        "type": "object",
        "properties": {
            "capacity_factor": {
                "type": "number",
                "description": (
                    "Capacity multiplier for the selected link. "
                    "A 50 percent increase means factor 1.5."
                ),
            },
        },
        "required": [
            "capacity_factor",
        ],
        "additionalProperties": False,
    },
}


EXPERIMENT_INSTRUCTIONS = """
You are Jeylan, the E-Routing network-planning assistant.

The user explicitly selected the Run Experiment operation.

Choose exactly one operation:

1. run_capacity_experiment:
   Use when the user explicitly identifies a physical edge ID
   to upgrade.

2. investigate_worst_bottleneck:
   Use when the user asks to find the highest-utilization
   physical link and test a capacity upgrade without
   providing an edge ID.

Rules:
- A 50 percent capacity increase means factor 1.5.
- The demand model defaults to uniform for explicit-edge experiments.
- Do not invent physical edge IDs.
- Do not invent congestion or utilization measurements.
- Do not request multiple experiments.
- Do not describe a heuristic selection as globally optimal.
- The selected network is supplied by the backend.
- The source analysis ID is supplied by the backend.
- Neither tool changes the active network workspace.
- If the request is outside the supported operations,
  do not invent another execution method.
"""


EXPLANATION_INSTRUCTIONS = """
You are Jeylan, the E-Routing network-planning assistant.

Explain a completed simulated network experiment using only
the supplied verified E-Routing evidence.

State:
- Baseline maximum congestion.
- Maximum congestion under fixed routing.
- Maximum congestion after routing recomputation.
- Relative changes compared with the baseline.

Lower congestion is better.

Do not assume that recomputing the routing scheme always
outperforms fixed routing.

If source_evidence is present:
- Identify the physical link selected from the baseline.
- State its source and target nodes when available.
- State its recorded baseline utilization.
- Identify the source analysis used for candidate selection.
- Explain that the selection criterion was highest baseline
  utilization, not global optimization of capacity investments.

These are modeled network-analysis results,
not measurements from a live production network.

Do not claim optimality, operational safety, or guaranteed
performance improvements.

Mention the persisted experiment ID.
Respond in the user's language.
"""


def _explicit_edge_ids(question: str) -> list[int]:
    """Extract explicitly named physical edge IDs."""

    return [
        int(value)
        for value in re.findall(
            r"\b(?:edge|link)\s*(?:id\s*)?#?\s*(\d+)\b",
            question,
            flags=re.IGNORECASE,
        )
    ]


def _parse_bottleneck_factor(arguments_json: str) -> float:
    """Validate a model-proposed bottleneck upgrade factor."""

    try:
        arguments = json.loads(arguments_json)

        if (
                not isinstance(arguments, dict)
                or set(arguments) != {"capacity_factor"}
        ):
            raise ValueError("Unexpected investigation arguments")

        # Reuse the same factor validation as explicit experiments.
        plan = CapacityExperimentPlan(
            edge_id=0,
            capacity_factor=arguments["capacity_factor"],
        )

        return plan.capacity_factor

    except (TypeError, ValueError) as exc:
        raise ValueError(
            "Jeylan proposed an invalid bottleneck investigation"
        ) from exc


def _fallback_explanation(
        measurements: dict,
        source_evidence: dict | None,
) -> str:
    """Provide measured results if the language-model explanation fails."""

    def fmt(value: object) -> str:
        if isinstance(value, (int, float)):
            return f"{value:,.3f}"
        return "Unavailable"

    baseline = measurements.get("baseline_max_congestion")
    fixed = measurements.get("fixed_routing_max_congestion")
    reoptimized = measurements.get("reoptimized_max_congestion")

    analysis_id = measurements["analysis_id"]
    edge_id = measurements.get("edge_id")
    factor = measurements.get("capacity_factor")

    lines = [
        "The simulated capacity experiment completed successfully.",
        "",
        f"Upgraded physical edge: {edge_id}",
        f"Capacity multiplier: {factor}",
        "",
        f"Baseline maximum congestion: {fmt(baseline)}",
        f"Fixed-routing maximum congestion: {fmt(fixed)}",
        f"Reoptimized maximum congestion: {fmt(reoptimized)}",
        "",
        f"Persisted analysis ID: {analysis_id}",
    ]

    if source_evidence is not None:
        candidate = source_evidence.get("candidate", {})

        lines.extend([
            "",
            "The link was selected using the highest baseline "
            "utilization criterion.",
            f"Selected baseline utilization: "
            f"{fmt(candidate.get('utilization'))}",
            f"Source analysis ID: "
            f"{source_evidence.get('source_analysis_id')}",
        ])

    return "\n".join(lines)


def answer_capacity_experiment(
        *,
        question: str,
        network_id: str,
        workspace: WorkspaceSnapshot,
        client: Any,
        model: str,
        source_analysis_id: str | None = None,
) -> dict:
    """
    Interpret and execute exactly one authorized capacity experiment.

    Supported operations:
    - Upgrade one explicitly specified physical edge.
    - Select the worst baseline bottleneck and upgrade it.

    Experiment execution is delegated to trusted E-Routing services.
    """

    if (
            not isinstance(question, str)
            or not 1 <= len(question.strip()) <= 2000
    ):
        raise ValueError(
            "Experiment request must contain 1–2000 characters"
        )

    # ---------------------------------------------------------
    # 1. Select exactly one model-facing operation.
    # ---------------------------------------------------------

    first = client.responses.create(
        model=model,
        instructions=EXPERIMENT_INSTRUCTIONS,
        input=[
            {
                "role": "user",
                "content": question,
            }
        ],
        tools=[
            CAPACITY_TOOL,
            BOTTLENECK_TOOL,
        ],
        tool_choice="required",
        parallel_tool_calls=False,
        store=False,
    )

    if first.status != "completed":
        raise RuntimeError(
            "Jeylan experiment planning failed"
        )

    calls = [
        item
        for item in first.output
        if item.type == "function_call"
    ]

    if len(calls) != 1:
        raise ValueError(
            "Expected exactly one experiment operation"
        )

    call = calls[0]
    source_evidence = None

    # ---------------------------------------------------------
    # 2A. Explicit physical-edge capacity experiment.
    # ---------------------------------------------------------

    if call.name == "run_capacity_experiment":
        try:
            plan = CapacityExperimentPlan.model_validate_json(
                call.arguments
            )
        except Exception as exc:
            raise ValueError(
                "Jeylan proposed an invalid capacity experiment"
            ) from exc

        mentioned_edges = _explicit_edge_ids(question)

        if (
                len(mentioned_edges) != 1
                or mentioned_edges[0] != plan.edge_id
        ):
            raise ValueError(
                "Specify exactly one physical edge ID for a "
                "direct capacity experiment, or request a "
                "worst-bottleneck investigation"
            )

        # Trusted backend execution.
        # No shell commands or network paths come from the model.
        result = run_capacity_experiment(
            network_id=network_id,
            plan=plan,
            workspace=workspace,
        )

        measurements = result.model_dump(
            mode="json"
        )

    # ---------------------------------------------------------
    # 2B. Evidence-driven bottleneck investigation.
    # ---------------------------------------------------------

    elif call.name == "investigate_worst_bottleneck":
        if source_analysis_id is None:
            raise ValueError(
                "Run a Capacity Planning analysis first to establish "
                "validated baseline link measurements"
            )

        factor = _parse_bottleneck_factor(
            call.arguments
        )

        # This service:
        # 1. Loads validated persisted baseline evidence.
        # 2. Selects the highest-utilization physical link.
        # 3. Executes exactly one bounded capacity experiment.
        investigation = run_bottleneck_investigation(
            network_id=network_id,
            source_analysis_id=source_analysis_id,
            workspace=workspace,
            capacity_factor=factor,
        )

        source_evidence = investigation["source_evidence"]
        measurements = investigation["experiment"]

    else:
        raise ValueError(
            "Unauthorized experiment operation"
        )

    # ---------------------------------------------------------
    # 3. Prepare evidence shared by both execution branches.
    # ---------------------------------------------------------

    analysis_id = measurements["analysis_id"]

    explanation_evidence = {
        "experiment": measurements,
        "source_evidence": source_evidence,
    }

    # ---------------------------------------------------------
    # 4. Generate an explanation without risking the results.
    # ---------------------------------------------------------

    try:
        second = client.responses.create(
            model=model,
            instructions=EXPLANATION_INSTRUCTIONS,
            input=[
                {
                    "role": "user",
                    "content": question,
                },
                {
                    "role": "user",
                    "content": (
                            "Verified E-Routing investigation evidence:\n"
                            + json.dumps(
                        explanation_evidence,
                        allow_nan=False,
                    )
                    ),
                },
            ],
            tools=[],
            store=False,
        )

        if (
                second.status != "completed"
                or not second.output_text.strip()
        ):
            answer = _fallback_explanation(
                measurements,
                source_evidence,
            )
        else:
            answer = second.output_text

    except Exception:
        # At this point the C++ experiment already completed.
        # Preserve its persisted results if the explanation request
        # encounters an OpenAI/network error.
        logger.exception(
            "Jeylan explanation failed after completed experiment %s",
            analysis_id,
        )

        answer = _fallback_explanation(
            measurements,
            source_evidence,
        )

    # ---------------------------------------------------------
    # 5. Return verified results and evidence provenance.
    # ---------------------------------------------------------

    evidence = [
        {
            "id": "E1",
            "source": "persisted_analysis_report",
            "scenario_id": analysis_id,
        }
    ]

    if source_evidence is not None:
        evidence.append({
            "id": "E2",
            "source": "persisted_analysis_report",
            "scenario_id": source_analysis_id,
        })

    return {
        "answer": answer,
        "network_id": network_id,
        "analysis_id": analysis_id,
        "analysis_type": "capacity_upgrade",
        "status": "completed",
        "experiment": measurements,
        "source_evidence": source_evidence,
        "evidence": evidence,
    }
