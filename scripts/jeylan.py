#!/usr/bin/env python3
"""Command-line interface for the Jeylan E-Routing assistant."""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

from e_routing.assistant.jeylan import JeylanAssistant
from e_routing.assistant.scenario_bridge import ScenarioBridge
from e_routing.scenarios.runtime import TrustedRun


def main() -> int:
    parser = argparse.ArgumentParser(description="Jeylan v0.1 — read-only E-Routing analyst")
    parser.add_argument("--results", type=Path, default=Path("results/canonical"),
                        help="Directory containing canonical JSON results")
    parser.add_argument("--schema", type=Path, default=Path("schema/analysis-result.schema.json"))
    parser.add_argument("--run", required=True, help="A single existing canonical run ID")
    parser.add_argument("--model", default=os.getenv("JEYLAN_MODEL", "gpt-5.4-mini"))
    parser.add_argument("--depth", choices=("executive", "operational", "technical"), default="operational")
    parser.add_argument("--ask", help="One-shot question; omit for an interactive session")
    parser.add_argument("--show-evidence", action="store_true", help="Print the tools used after each answer")
    parser.add_argument("--engine", type=Path)
    parser.add_argument("--graph", type=Path)
    parser.add_argument("--routing-table", type=Path)
    parser.add_argument("--scenario-solver")
    parser.add_argument("--scenario-workspace", type=Path, default=Path("results/scenarios"))
    parser.add_argument("--scenario-seed", type=int, default=42)

    args = parser.parse_args()

    if not os.getenv("OPENAI_API_KEY"):
        parser.error("Set OPENAI_API_KEY in your environment before running Jeylan")
    try:
        from openai import OpenAI

        from e_routing.analysis.queries import AnalysisAPI
    except ImportError as exc:
        parser.error(f"Failed to load E-Routing analysis dependencies: {exc}")

    import re

    if not re.fullmatch(r"[A-Za-z0-9_-]+", args.run):
        parser.error("Invalid run ID")

    scenario_path = args.results / f"{args.run}.json"

    try:
        with scenario_path.open(encoding="utf-8") as f:
            scenario_data = json.load(f)
    except (OSError, ValueError) as exc:
        parser.error(f"Unable to load scenario data: {exc}")

    if scenario_data.get("analysis_metadata", {}).get("run_id") != args.run:
        parser.error("Canonical result run ID does not match selected run")

    scenario_bridge = None
    scenario_bridge = None

    scenario_options = [
        args.engine,
        args.graph,
        args.routing_table,
        args.scenario_solver,
    ]

    if any(option is not None for option in scenario_options):
        if not all(option is not None for option in scenario_options):
            parser.error("--engine, --graph, --routing-table and --scenario-solver must be provided together")

        scenario_bridge = ScenarioBridge(
            TrustedRun(
                executable=args.engine,
                graph=args.graph,
                routing_table=args.routing_table,
                solver_label=args.scenario_solver,
                workspace=args.scenario_workspace,
                seed=args.scenario_seed,
            )
        )

    agent = JeylanAssistant(
        AnalysisAPI(args.results, args.schema),
        OpenAI(),
        args.run,
        args.model,
        args.depth,
        scenario_data=scenario_data,
        scenario_bridge=scenario_bridge
    )

    def answer(question: str) -> None:
        print(f"\nJeylan: {agent.ask(question)}")
        if args.show_evidence:
            for e in agent.last_evidence:
                if e.tool == "evaluate_persisted_scenario":
                    print(f"  [{e.id}] {e.tool}({e.arguments}) "
                          f"scenario={e.data.get('scenario_id')}")
                else:
                    print(f"  [{e.id}] {e.tool}({e.arguments}) "
                          f"run={args.run}")

    try:
        if args.ask:
            answer(args.ask)
            return 0
        print(f"Jeylan v0.1 | run={args.run} | depth={args.depth} | /reset, /exit")
        while True:
            try:
                question = input("\nYou: ").strip()
            except EOFError:
                break
            if question.lower() in {"/exit", "/quit"}:
                break
            if question.lower() == "/reset":
                agent.reset()
                print("Conversation cleared.")
                continue
            if question:
                answer(question)
    except (RuntimeError, ValueError, OSError) as exc:
        print(f"Jeylan error: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
