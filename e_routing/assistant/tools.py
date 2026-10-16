
"""Validated deterministic tools available to Jeylan."""

from __future__ import annotations

from typing import Any

from e_routing.analysis.scenarios import assess_link_maintenance

TOOL_DEFINITIONS: list[dict[str, Any]] = [
    {
        "type": "function", "name": "get_analysis_summary", "strict": True,
        "description": "Retrieve the authoritative analysis summary for the currently selected E-Routing run.",
        "parameters": {"type": "object", "properties": {}, "required": [], "additionalProperties": False},
    },
    {
        "type": "function", "name": "get_failure_impact", "strict": True,
        "description": "Retrieve recorded failure impact for one physical edge ID. No new simulation is performed.",
        "parameters": {
            "type": "object",
            "properties": {"failed_edge_id": {"type": "integer",
                                              "description": "Physical edge identifier from analysis evidence."}},
            "required": ["failed_edge_id"], "additionalProperties": False,
        },
    },
    {
        "type": "function", "name": "rank_critical_links", "strict": True,
        "description": "Rank network links using the available deterministic congestion-increase metric.",
        "parameters": {
            "type": "object",
            "properties": {
                "metric": {"type": "string", "enum": ["congestion_increase_factor"]},
                "limit": {"type": "integer", "description": "Maximum number of ranked links (1 through 20)."},
            },
            "required": ["metric", "limit"], "additionalProperties": False,
        },
    },

    {
        "type": "function",
        "name": "compare_link_failures",
        "strict": True,
        "description": (
            "Compare already evaluated, independent N-1 link failures "
            "for 2 to 10 physical edges. This is not a simultaneous "
            "multi-link failure simulation."
        ),
        "parameters": {
            "type": "object",
            "properties": {
                "edge_ids": {
                    "type": "array",
                    "items": {"type": "integer"},
                    "description": "Physical edge IDs to compare."
                }
            },
            "required": ["edge_ids"],
            "additionalProperties": False
        }
    },

    {
        "type": "function",
        "name": "assess_link_maintenance",
        "strict": True,
        "description": (
            "Assess planned maintenance of one physical link using "
            "recorded N-1 recovery evidence. Returns findings, "
            "concerns, and required operational checks. "
            "Does not authorize maintenance."
        ),
        "parameters": {
            "type": "object",
            "properties": {
                "edge_id": {
                    "type": "integer",
                    "description": "Physical edge ID proposed for maintenance."
                }
            },
            "required": ["edge_id"],
            "additionalProperties": False
        }
    },

    {
        "type": "function",
        "name": "evaluate_persisted_scenario",
        "strict": True,
        "description": (
            "Run a new baseline demand evaluation using a prevalidated, "
            "persisted E-Routing scheme. Optionally compute static Layer-1 "
            "single-link failure exposure. This does not perform routing "
            "optimization, link removal, recovery, or demand scaling."
        ),
        "parameters": {
            "type": "object",
            "properties": {
                "demand_model": {
                    "type": "string",
                    "enum": ["uniform", "gravity", "bimodal", "gaussian"]
                },
                "include_layer1": {
                    "type": "boolean"
                },

                "demand_scale": {
                    "type": "number",
                    "description": "Demand multiplier, e.g. 1.3 for 30 percent growth"
                },
                "capacity_upgrades": {
                    "type": "array",
                    "description": "Hypothetical capacity upgrades to canonical physical link IDs",
                    "items": {
                        "type": "object",
                        "properties": {
                            "edge_id": {"type": "integer"},
                            "capacity_factor": {"type": "number"}
                        },
                        "required": ["edge_id", "capacity_factor"],
                        "additionalProperties": False
                    }
                }
            },
            "required": ["demand_model", "include_layer1", "demand_scale", "capacity_upgrades"],
            "additionalProperties": False
        }
    },

]

def require_exact_keys(args: Any, required: set[str]) -> dict[str, Any]:
    if not isinstance(args, dict) or set(args) != required:
        raise ValueError(f"Expected arguments: {', '.join(sorted(required)) or '(none)'}")
    return args


class ToolDispatcher:
    def __init__(self, analysis_api: Any, run_id: str, scenario_data: dict[str, Any] | None = None, scenario_bridge: Any | None = None):
        self.api = analysis_api
        self.run_id = run_id
        self.scenario_data = scenario_data
        self.scenario_bridge = scenario_bridge

    def dispatch(self, name: str, args: dict[str, Any]) -> Any:
        if name == "get_analysis_summary":
            require_exact_keys(args, set())
            return self.api.get_analysis_summary(self.run_id)

        if name == "get_failure_impact":
            require_exact_keys(args, {"failed_edge_id"})
            edge_id = args["failed_edge_id"]
            if type(edge_id) is not int or edge_id < 0:
                raise ValueError("failed_edge_id must be a nonnegative integer")
            return self.api.get_failure_impact(self.run_id, edge_id)

        if name == "rank_critical_links":
            require_exact_keys(args, {"metric", "limit"})
            if args["metric"] != "congestion_increase_factor":
                raise ValueError("Unsupported ranking metric")
            limit = args["limit"]
            if type(limit) is not int or not 1 <= limit <= 20:
                raise ValueError("limit must be between 1 and 20")
            return self.api.rank_critical_links(self.run_id, args["metric"], limit)

        if name == "compare_link_failures":
            require_exact_keys(args, {"edge_ids"})
            edge_ids = args["edge_ids"]
            if not isinstance(edge_ids, list) or not 2 <= len(edge_ids) <= 10:
                raise ValueError("Provide between 2 and 10 edge IDs")
            if any(type(edge) is not int or edge < 0 for edge in edge_ids):
                raise ValueError("All edge IDs must be nonnegative integers")
            if len(set(edge_ids)) != len(edge_ids):
                raise ValueError("Duplicate edge IDs are not allowed")

            return {
                "scenario_type": "independent_n1_failure_comparison",
                "run_id": self.run_id,
                "results": [
                    {"edge_id": edge_id, "impact": self.api.get_failure_impact(self.run_id, edge_id)}
                    for edge_id in edge_ids
                ],
                "limitations": (
                    "Each edge was evaluated as an independent single-link failure. "
                    "These results do not represent simultaneous failures or failure probabilities."
                ),
            }

        if name == "assess_link_maintenance":
            require_exact_keys(args, {"edge_id"})
            edge_id = args["edge_id"]
            if type(edge_id) is not int or edge_id < 0:
                raise ValueError("edge_id must be a nonnegative integer")
            if self.scenario_data is None:
                raise ValueError("Canonical scenario data unavailable")
            return assess_link_maintenance(self.scenario_data, edge_id)

        if name == "evaluate_persisted_scenario":
            require_exact_keys(args, {"demand_model", "include_layer1", "demand_scale", "capacity_upgrades"})
            if self.scenario_bridge is None:
                raise ValueError("Persisted-scheme scenario execution is not configured")
            return self.scenario_bridge.evaluate(args)

        raise ValueError(f"Unknown/unauthorized tool: {name}")
