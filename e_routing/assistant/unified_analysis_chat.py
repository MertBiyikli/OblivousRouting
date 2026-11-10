
"""One authorized, evidence-grounded Jeylan network analysis."""

from __future__ import annotations

import json
import logging
import re
from pathlib import Path
from typing import Any

from fastapi.concurrency import run_in_threadpool
from pydantic import TypeAdapter, ValidationError

from e_routing.api.workspace_models import WorkspaceSnapshot
from e_routing.assistant.analysis_orchestrator import (
    run_analysis_plan,
)
from e_routing.assistant.analysis_plan import AnalysisOperation

from e_routing.assistant.composite_plans import (
    CompositeAnalysisPlan,
    CompositeCapacityStep,
    CompositeLinkFailureStep,
)
from e_routing.assistant.composite_executor import (
    execute_composite_analysis,
)


logger = logging.getLogger(__name__)

PLAN_ADAPTER = TypeAdapter(AnalysisOperation)

DEMAND_MODELS = [
    "uniform", "gravity", "bimodal", "gaussian"
]


def _tool(name: str, description: str, properties: dict) -> dict:
    return {
        "type": "function",
        "name": name,
        "description": description,
        "strict": True,
        "parameters": {
            "type": "object",
            "properties": properties,
            "required": list(properties),
            "additionalProperties": False,
        },
    }


DEMAND_PROPERTY = {
    "type": "string",
    "enum": DEMAND_MODELS,
}



def _composite_step_schema(
        operation: str,
        properties: dict,
) -> dict:
    """Strict schema for one supported composite step."""

    fields = {
        "operation": {
            "type": "string",
            "enum": [operation],
        },
        **properties,
    }

    return {
        "type": "object",
        "properties": fields,
        "required": list(fields),
        "additionalProperties": False,
    }


COMPOSITE_STEP_SCHEMA = {
    "anyOf": [
        _composite_step_schema(
            "traffic_growth",
            {
                "traffic_multiplier": {
                    "type": "number",
                    "exclusiveMinimum": 0,
                    "maximum": 1000,
                },
            },
        ),
        _composite_step_schema(
            "capacity_upgrade",
            {
                "upgrades": {
                    "type": "array",
                    "minItems": 1,
                    "maxItems": 1,
                    "items": {
                        "type": "object",
                        "properties": {
                            "edge_id": {
                                "type": "integer",
                            },
                            "capacity_factor": {
                                "type": "number",
                            },
                        },
                        "required": [
                            "edge_id",
                            "capacity_factor",
                        ],
                        "additionalProperties": False,
                    },
                },
            },
        ),
        _composite_step_schema(
            "n1_resilience",
            {
                "execution_mode": {
                    "type": "string",
                    "enum": ["fast", "standard"],
                },
            },
        ),
        _composite_step_schema(
            "link_failure",
            {
                "failed_edge_id": {
                    "type": "integer",
                },
            },
        ),
    ],
}


TOOLS = [
    _tool(
        "capacity_upgrade",
        "Simulate one explicitly specified physical-link capacity "
        "upgrade and compare baseline, fixed and reoptimized routing.",
        {
            "demand_model": DEMAND_PROPERTY,
            "upgrades": {
                "type": "array",
                "items": {
                    "type": "object",
                    "properties": {
                        "edge_id": {"type": "integer"},
                        "capacity_factor": {"type": "number"},
                    },
                    "required": ["edge_id", "capacity_factor"],
                    "additionalProperties": False,
                },
            },
        },
    ),
    _tool(
        "link_failure",
        "Simulate the failure of one physical link explicitly "
        "identified by the user.",
        {
            "demand_model": DEMAND_PROPERTY,
            "failed_edge_id": {"type": "integer"},
        },
    ),
    _tool(
        "n1_resilience",
        "Evaluate all single-physical-link failures using "
        "the current network workspace.",
        {
            "demand_model": DEMAND_PROPERTY,
            "execution_mode": {
                "type": "string",
                "enum": ["fast", "standard"],
            },
            "seed": {"type": "integer"},
        },
    ),
    _tool(
        "traffic_growth",
        "Evaluate one relative traffic-demand multiplier "
        "against the current workspace using fixed routing.",
        {
            "demand_model": DEMAND_PROPERTY,
            "traffic_multiplier": {"type": "number"},
        },
    ),

    _tool(
        "composite_analysis",
        "Execute a bounded, sequential network investigation. "
        "Supports Traffic Growth, one explicitly specified "
        "Capacity Upgrade, N-1 Resilience and Link Failure. "
        "Use when the user requests two or more analysis steps. "
        "Preserve the requested order. Never invent edge IDs "
        "or capacity candidates.",
        {
            "demand_model": DEMAND_PROPERTY,
            "seed": {
                "type": "integer",
                "enum": [42],
            },
            "steps": {
                "type": "array",
                "minItems": 2,
                "maxItems": 6,
                "items": COMPOSITE_STEP_SCHEMA,
            },
        },
    ),
]



INSTRUCTIONS = """
You are Jeylan, the E-Routing network-analysis assistant.

The user explicitly requested a simulated network analysis.

Select exactly one execution tool.

For one operation, use one of:
- capacity_upgrade
- link_failure
- n1_resilience
- traffic_growth

For two or more sequential operations, use:
- composite_analysis

Composite analysis supports:
- Traffic Growth
- Capacity Upgrade
- N-1 Resilience
- Link Failure

Preserve the order and purpose of the user's request.

The network and existing workspace are provided by the
backend. Never invent or override them.

A 30 percent traffic increase means multiplier 1.3.
A 50 percent capacity increase means factor 1.5.

Traffic Growth and Capacity Upgrade change only the
simulated working workspace for subsequent steps.

N-1 Resilience and Link Failure evaluate the current
simulated configuration without permanently changing it.

Capacity upgrades and individual link failures require
explicit physical edge IDs supplied by the user.
Never invent edge IDs or capacity upgrade candidates.

The current composite implementation supports:
- 2 to 6 steps
- at most one Traffic Growth step
- at most one Capacity Upgrade step
- at most one edge per Capacity Upgrade step
- supported demand models only
- FAST or STANDARD N-1 execution
- seed 42

Use uniform demand, FAST N-1 mode and seed 42 unless
another supported option is explicitly requested.
Do not silently substitute an unsupported seed.

Do not invent traffic multipliers, capacity factors,
measurements or execution results.

Do not propose autonomous optimization or candidate
search as though it were an authorized analysis.

If a request requires unsupported operations or
ambiguous physical edge IDs, do not call a tool.

These are simulations. Do not claim that production
routing or the active Scenario Studio workspace changes.

A composite execution may return partial results.
Never describe a partial result as fully completed.
"""



EXPLANATION_INSTRUCTIONS = """
Explain the completed E-Routing analysis using only the
supplied validated measurements.

Distinguish fixed routing from recomputed routing.
For N-1 analysis, distinguish disconnected failures from
successful rerouting. Never imply that all failures succeeded
when disconnected outcomes exist.

Lower maximum congestion is better.
Congestion is a load/capacity ratio.

Traffic Growth uses fixed routing and an effective demand
scale derived from the existing workspace.

Do not claim optimality or operational safety.
These are simulated results, not live telemetry.

Mention the persisted analysis ID.
Respond in the user's language.

Present congestion values to two decimal places.
Present relative changes as percentages with two decimal places.
Present demand volumes with at most two decimal places.

Do not expose floating-point representation artifacts.
Do not round the persisted report or change the underlying
measurements.

Distinguish the requested relative traffic multiplier
from the effective demand scale of the workspace.
"""


COMPOSITE_EXPLANATION_INSTRUCTIONS = """
You are explaining a completed or partially completed
E-Routing composite network investigation.

Use only the supplied verified measurements.

Explain the operations in their actual execution order.
Mention the persisted analysis IDs.

Clearly distinguish:
- Traffic Growth
- Capacity Upgrade
- N-1 Resilience
- Link Failure

A lower maximum congestion value is better.

For capacity comparisons, distinguish:
- baseline congestion changes
- post-failure congestion changes
- connectivity improvements
- connectivity regressions
- execution failures

Do not treat disconnected failures as numerical congestion.
Do not claim optimality or production safety.

If a comparison contains truncated edge details, explain
only the supplied individual edges. You may report complete
aggregate counts when supplied.

If the overall status is partial, clearly identify what
completed and what did not.

When comparison_status is failed, do not claim that a
cross-analysis comparison was verified.

Never invent unreported measurements or edge rankings.
Mention that the analysis is simulated, not live telemetry.

Respond in the user's language.

Display congestion values with two decimal places.
Display relative changes as percentages with two decimals.
Avoid floating-point representation artifacts.
"""



def _explicit_edge_ids(question: str) -> list[int]:
    return [
        int(value)
        for value in re.findall(
            r"\b(?:edge|link|kante)\s*"
            r"(?:id\s*)?#?\s*(\d+)\b",
            question,
            flags=re.IGNORECASE,
        )
    ]


def _reject_known_composite_request(question: str) -> None:
    """Conservatively block the initial traffic-plus-resilience case."""
    traffic = re.search(
        r"\b(?:traffic|demand|verkehr|nachfrage)\b",
        question,
        re.IGNORECASE,
    )
    resilience = re.search(
        r"(?:n[-−]1|single[- ]link|resilien|"
        r"failure|ausfall)",
        question,
        re.IGNORECASE,
    )

    if traffic and resilience:
        raise ValueError(
            "Combined Traffic Growth and resilience investigations "
            "are not supported yet. Run one analysis at a time."
        )


def _validated_plan(call: Any, question: str):
    allowed = {
        "capacity_upgrade",
        "link_failure",
        "n1_resilience",
        "traffic_growth",
    }

    if call.name not in allowed:
        raise ValueError("Unauthorized analysis operation")

    try:
        arguments = json.loads(call.arguments)
        if not isinstance(arguments, dict):
            raise ValueError("Tool arguments must be an object")

        plan = PLAN_ADAPTER.validate_python({
            "kind": call.name,
            **arguments,
        })

    except (ValueError, TypeError, ValidationError) as exc:
        raise ValueError(
            "Jeylan proposed an invalid analysis operation"
        ) from exc

    if plan.kind in {"capacity_upgrade", "link_failure"}:
        mentioned = _explicit_edge_ids(question)

        expected = (
            plan.upgrades[0].edge_id
            if plan.kind == "capacity_upgrade"
            else plan.failed_edge_id
        )

        if len(mentioned) != 1 or mentioned[0] != expected:
            raise ValueError(
                "Specify exactly one physical edge ID for this "
                "analysis. Jeylan must not invent edge IDs."
            )

    return plan


def _validated_composite_plan(
        call: Any,
        question: str,
) -> CompositeAnalysisPlan:
    """Validate a model-proposed composite execution plan."""

    if call.name != "composite_analysis":
        raise ValueError("Unauthorized composite operation")

    try:
        arguments = json.loads(call.arguments)

        if not isinstance(arguments, dict):
            raise ValueError("Tool arguments must be an object")

        if "kind" in arguments:
            raise ValueError("Model cannot override plan kind")

        plan = CompositeAnalysisPlan.model_validate({
            "kind": "composite",
            **arguments,
        })

    except (ValueError, TypeError, ValidationError) as exc:
        raise ValueError(
            "Jeylan proposed an invalid composite plan"
        ) from exc

    # Only user-explicit physical edges may be modified or failed.
    mentioned = set(_explicit_edge_ids(question))
    requested = set()

    for step in plan.steps:
        if isinstance(step, CompositeCapacityStep):
            requested.update(
                item.edge_id for item in step.upgrades
            )

        elif isinstance(step, CompositeLinkFailureStep):
            requested.add(step.failed_edge_id)

    if not requested.issubset(mentioned):
        raise ValueError(
            "All capacity upgrades and individual link failures "
            "require physical edge IDs explicitly named by "
            "the user. Jeylan must not invent edge IDs."
        )

    return plan


def _fallback_answer(result: dict) -> str:
    return (
            "The simulated network analysis completed successfully.\n\n"
            f"Analysis type: {result['analysis_type']}\n"
            f"Persisted analysis ID: {result['analysis_id']}\n\n"
            "Verified summary:\n"
            + json.dumps(
        result["summary"],
        ensure_ascii=False,
        indent=2,
        allow_nan=False,
    )
    )


def _composite_evidence(result: dict) -> dict:
    """Build bounded, verified evidence for Jeylan."""

    steps = []

    for step in result["completed_steps"]:
        steps.append({
            "step_index": step["step_index"],
            "operation": step["operation"],
            "analysis_id": step["analysis_id"],
            "plan": step["plan"],
            "summary": step["summary"],
            "input_workspace": step["input_workspace"],
            "output_workspace": step["output_workspace"],
        })

    comparisons = []

    for comparison in result.get("comparisons", [])[:5]:
        rows = comparison.get("edge_comparisons", [])

        # Retain aggregate counts, but bound detailed rows.
        comparisons.append({
            "comparison_type": comparison["comparison_type"],
            "before_analysis_id":
                comparison["before_analysis_id"],
            "after_analysis_id":
                comparison["after_analysis_id"],
            "before_baseline_max_congestion":
                comparison["before_baseline_max_congestion"],
            "after_baseline_max_congestion":
                comparison["after_baseline_max_congestion"],
            "relative_baseline_congestion_change":
                comparison["relative_baseline_congestion_change"],
            "counts": comparison["counts"],
            "edge_comparisons": rows[:10],
            "edge_comparisons_total": len(rows),
            "edge_comparisons_truncated": len(rows) > 10,
        })

    return {
        "network_id": result["network_id"],
        "status": result["status"],
        "kind": "composite",
        "plan": result["plan"],
        "completed_steps": steps,
        "comparisons": comparisons,
        "comparison_status": result.get(
            "comparison_status", "not_available"
        ),
        "failed_step_index": result.get("failed_step_index"),
        "failed_operation": result.get("failed_operation"),
        "error": result.get("error"),
        "workspace_modified": False,
    }


def _fallback_composite_answer(result: dict) -> str:
    """Deterministic explanation when model explanation fails."""

    lines = [
        "Simulated composite network investigation",
        f"Status: {result['status']}",
        "",
        "Persisted analyses:",
    ]

    for step in result["completed_steps"]:
        lines.append(
            f"- Step {step['step_index'] + 1}: "
            f"{step['operation']} "
            f"({step['analysis_id']})"
        )

    if result.get("comparisons"):
        lines.append("")
        lines.append(
            f"Verified comparisons: {len(result['comparisons'])}"
        )

        for comparison in result["comparisons"]:
            lines.append(
                json.dumps({
                    "before_analysis_id":
                        comparison["before_analysis_id"],
                    "after_analysis_id":
                        comparison["after_analysis_id"],
                    "counts": comparison["counts"],
                }, ensure_ascii=False, allow_nan=False)
            )

    if result["status"] != "completed":
        lines.append("")
        lines.append(
            "The investigation did not fully complete. "
            "Only the listed analyses are confirmed completed."
        )

        if result.get("error"):
            lines.append(result["error"])

    return "\n".join(lines)


async def answer_unified_analysis(
        *,
        question: str,
        network_id: str,
        workspace: WorkspaceSnapshot,
        client: Any,
        model: str,
        reports_root: Path | str = "results/scenarios",
) -> dict:
    """Plan and execute at most one authorized network analysis."""

    if not isinstance(question, str) or not 1 <= len(
            question.strip()
    ) <= 2000:
        raise ValueError("Analysis request must contain 1-2000 characters")


    first = await run_in_threadpool(
        client.responses.create,
        model=model,
        instructions=INSTRUCTIONS,
        input=[{"role": "user", "content": question}],
        tools=TOOLS,
        tool_choice="auto",
        parallel_tool_calls=False,
        store=False,
    )

    if first.status != "completed":
        raise RuntimeError("Jeylan analysis planning failed")

    calls = [
        item for item in first.output
        if item.type == "function_call"
    ]

    if len(calls) != 1:
        raise ValueError("Request exactly one supported execution tool")


    call = calls[0]

    # Prevent incomplete single-operation execution when
    # the user explicitly requests Traffic Growth + resilience.
    if call.name != "composite_analysis":
        traffic_requested = re.search(
            r"\b(?:traffic|demand|verkehr|nachfrage)\b",
            question,
            re.IGNORECASE,
        )

        resilience_requested = re.search(
            r"(?:n[-−]1|single[- ]link|resilien|"
            r"failure|ausfall)",
            question,
            re.IGNORECASE,
        )

        if traffic_requested and resilience_requested:
            raise ValueError(
                "The request contains multiple analyses. "
                "Jeylan must select a composite_analysis tool "
                "instead of executing only one operation."
            )

    if call.name == "composite_analysis":
        plan = _validated_composite_plan(call, question)
    else:
        plan = _validated_plan(call, question)



    if isinstance(plan, CompositeAnalysisPlan):
        # Reuse the validated F3-A/F3-B execution engine.
        result = await execute_composite_analysis(
            network_id=network_id,
            plan=plan,
            workspace=workspace,
            reports_root=reports_root,
        )

        if result.get("status") not in {
            "completed", "partial"
        }:
            raise RuntimeError(
                "Composite analysis returned invalid status"
            )

        if result.get("network_id") != network_id:
            raise RuntimeError(
                "Composite analysis network identity mismatch"
            )

        evidence = _composite_evidence(result)

        try:
            second = await run_in_threadpool(
                client.responses.create,
                model=model,
                instructions=COMPOSITE_EXPLANATION_INSTRUCTIONS,
                input=[
                    {
                        "role": "user",
                        "content": question,
                    },
                    {
                        "role": "user",
                        "content": (
                                "Verified E-Routing composite evidence:\n"
                                + json.dumps(
                            evidence,
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
                answer = _fallback_composite_answer(result)
            else:
                answer = second.output_text

        except Exception:
            logger.exception(
                "Jeylan composite explanation failed"
            )
            answer = _fallback_composite_answer(result)

        persisted_evidence = [
            {
                "id": f"E{index + 1}",
                "source": "persisted_analysis_report",
                "scenario_id": step["analysis_id"],
            }
            for index, step in enumerate(
                result["completed_steps"]
            )
        ]

        return {
            "answer": answer,
            "network_id": network_id,
            "analysis_id": None,
            "analysis_type": "composite",
            "status": result["status"],
            "plan": result["plan"],
            "completed_steps": result["completed_steps"],
            "comparisons": result.get("comparisons", []),
            "comparison_status": result.get(
                "comparison_status", "not_available"
            ),
            "failed_step_index": result.get(
                "failed_step_index"
            ),
            "failed_operation": result.get(
                "failed_operation"
            ),
            "error": result.get("error"),
            "effective_workspace":
                result["effective_workspace"],
            "workspace_modified": False,
            "evidence": persisted_evidence,
        }


    # The trusted orchestrator performs at most one execution.
    result = await run_analysis_plan(
        network_id=network_id,
        plan=plan,
        workspace=workspace,
        reports_root=reports_root,
    )

    if result["status"] != "completed":
        raise RuntimeError("Analysis did not complete")

    # Do not send the potentially large topology/report payload
    # back to the language model. Only send bounded measurements.
    evidence = {
        "network_id": result["network_id"],
        "analysis_id": result["analysis_id"],
        "analysis_type": result["analysis_type"],
        "plan": result["plan"],
        "summary": result["summary"],
    }

    try:
        second = await run_in_threadpool(
            client.responses.create,
            model=model,
            instructions=EXPLANATION_INSTRUCTIONS,
            input=[
                {"role": "user", "content": question},
                {
                    "role": "user",
                    "content": (
                            "Verified E-Routing results:\n"
                            + json.dumps(evidence, allow_nan=False)
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
            answer = _fallback_answer(result)
        else:
            answer = second.output_text

    except Exception:
        # Execution already completed; never discard its result.
        logger.exception(
            "Jeylan explanation failed after analysis %s",
            result["analysis_id"],
        )
        answer = _fallback_answer(result)

    return {
        "answer": answer,
        "network_id": result["network_id"],
        "analysis_id": result["analysis_id"],
        "analysis_type": result["analysis_type"],
        "status": "completed",
        "summary": result["summary"],
        "plan": result["plan"],
        "evidence": [{
            "id": "E1",
            "source": "persisted_analysis_report",
            "scenario_id": result["analysis_id"],
        }],
    }
