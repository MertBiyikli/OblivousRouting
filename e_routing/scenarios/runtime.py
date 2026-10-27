#!/usr/bin/env python3
"""Deterministic E-Routing persisted-scheme scenario runner (v0.1).

Only trusted application code supplies the binary, graph and routing-file paths.
Model-visible requests may select a demand model and Layer-1 analysis, nothing else.
"""
from __future__ import annotations

import argparse
import json
import math
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Any
from dataclasses import replace
from uuid import uuid4

from e_routing.scenarios.fixed_routing import transform_fixed_routing

ALLOWED_MODELS = frozenset({'uniform', 'gravity', 'bimodal', 'gaussian'})



def execute_comparison(
        request: ScenarioRequest,
        trusted: TrustedRun,
) -> dict[str, Any]:
    if request.scenario_type != "capacity_planning":
        raise ValueError("Comparison requires capacity_planning")

    if request.include_layer1:
        raise ValueError("Fast comparison does not support Layer-1 analysis")

    if request.demand_model not in {"uniform", "bimodal", "gaussian"}:
        raise ValueError("Demand model is not eligible for fast comparison")

    baseline_request = replace(
        request,
        scenario_id="baseline-" + uuid4().hex,
        scenario_type="baseline",
        demand_scale=1.0,
        capacity_upgrades=(),
        include_layer1=False,
    )

    # The only invocation of the C++ solver.
    baseline = execute(baseline_request, trusted)

    reference = transform_fixed_routing(
        baseline,
        scenario_id=request.scenario_id + "-reference",
        demand_scale=request.demand_scale,
    )

    final = transform_fixed_routing(
        baseline,
        scenario_id=request.scenario_id,
        demand_scale=request.demand_scale,
        capacity_upgrades=request.capacity_upgrades,
    )

    final["comparison"] = {
        "baseline": baseline,
        "reference": reference,
    }

    return final


@dataclass(frozen=True)
class ScenarioRequest:
    scenario_id: str
    scenario_type: str
    demand_model: str
    include_layer1: bool
    demand_scale: float = 1.0
    capacity_upgrades: tuple[tuple[int, float], ...] = ()

    @classmethod
    def from_dict(cls, raw: dict[str, Any]) -> ScenarioRequest:
        if not isinstance(raw, dict):
            raise ValueError('Scenario request must be an object')
        required = {'scenario_id', 'scenario_type', 'demand_model', 'include_layer1'}
        allowed = required | {'demand_scale', 'capacity_upgrades'}
        if not required <= set(raw) or not set(raw) <= allowed:
            raise ValueError(
                f'Missing required fields: {sorted(required - set(raw))}; '
                f'unknown fields: {sorted(set(raw) - allowed)}'
            )
        sid = raw['scenario_id']
        if not isinstance(sid, str) or not 1 <= len(sid) <= 80 or not all(c.isascii() and (c.isalnum() or c in "-_") for c in sid):
            raise ValueError('scenario_id must be 1..80 alphanumeric, dash or underscore characters')
        if raw['scenario_type'] not in {'baseline', 'capacity_planning'}:
            raise ValueError('Only baseline and capacity_planning scenarios are supported')
        if raw['demand_model'] not in ALLOWED_MODELS:
            raise ValueError('Unsupported demand model')
        if type(raw['include_layer1']) is not bool:
            raise ValueError('include_layer1 must be boolean')
        scale = raw.get('demand_scale', 1.0)
        if isinstance(scale, bool) or not isinstance(scale, (float, int)) or not math.isfinite(scale) or not 0 < scale <= 100:
            raise ValueError('demand_scale must be finite and between 0 and 100')
        upgrades = raw.get('capacity_upgrades', [])
        if not isinstance(upgrades, list) or len(upgrades) > 100:
            raise ValueError('capacity_upgrades must be a list of at most 100 upgrades')
        parsed = []
        seen = set()
        for upgrade in upgrades:
            if not isinstance(upgrade, dict) or set(upgrade) != {'edge_id', 'capacity_factor'}:
                raise ValueError('Capacity upgrade requires edge_id and capacity_factor')
            edge, factor = upgrade['edge_id'], upgrade['capacity_factor']
            if (type(edge) is not int
                    or edge < 0
                    or edge in seen):
                raise ValueError('Capacity edge IDs must be unique nonnegative integers')
            if isinstance(factor, bool) or not isinstance(factor, (float, int)) or not math.isfinite(factor) or not 1 <= factor <= 100:
                raise ValueError('capacity_factor must be finite and between 1 and 100')
            seen.add(edge)
            parsed.append((edge, float(factor)))
        if raw['scenario_type'] == 'baseline' and (float(scale) != 1.0 or parsed):
            raise ValueError('Use capacity_planning for changed demand or capacity')
        return cls(sid, raw['scenario_type'], raw['demand_model'], raw['include_layer1'], float(scale), tuple(parsed))


@dataclass(frozen=True)
class TrustedRun:
    executable: Path
    graph: Path
    routing_table: Path
    solver_label: str
    workspace: Path
    seed: int = 42


def execute(request: ScenarioRequest, trusted: TrustedRun, timeout_seconds: int = 120) -> dict[str, Any]:
    if not isinstance(request, ScenarioRequest):
        raise TypeError('Expected validated ScenarioRequest')

    # TODO: put the solver labels in a config file or environment variable instead of hardcoding them here.
    if trusted.solver_label not in {'electrical', 'electrical_naive', 'raecke_frt', 'raecke_ckr', 'raecke_mst', 'cohen', 'expander', 'expander_mwu', 'semi_electrical', 'semi_tree', 'semi_expander'}:
        raise ValueError('Invalid trusted solver label')

    # Resolve trusted paths and prepare the command line for the trusted binary.
    binary = trusted.executable.resolve(strict=True)
    graph = trusted.graph.resolve(strict=True)
    table = trusted.routing_table.resolve(strict=True)
    workspace = trusted.workspace.resolve()
    workspace.mkdir(parents=True, exist_ok=True)
    result_path = workspace / (request.scenario_id + '.json')


    # Use an argv array and a trusted absolute executable path: never shell=True.
    command = [str(binary), 'solve', '--solver', trusted.solver_label,
               '--graph', str(graph),
               '--load-routing-table', str(table),
               '--demand', request.demand_model,
               '--seed', str(trusted.seed),
               '--output', str(result_path)]

    command.extend(['--demand-scale', str(request.demand_scale)])

    # Add capacity upgrades and Layer-1 analysis flags if requested.
    for edge, factor in request.capacity_upgrades:
        command.extend(['--capacity-upgrade', f'{edge}:{factor}'])
    if request.include_layer1:
        command.append('--failure-analysis')

    # Run the trusted binary and check for errors, timeouts, and result shape.
    try:
        completed = subprocess.run(command, text=True, capture_output=True, timeout=timeout_seconds, check=False)
    except subprocess.TimeoutExpired as exc:
        raise RuntimeError('E-Routing scenario evaluation timed out') from exc
    if completed.returncode != 0:
        raise RuntimeError(f'E-Routing exited with {completed.returncode}: {completed.stderr[-1500:]}')
    data = json.loads(result_path.read_text(encoding='utf-8'))
    results = data.get('solver_results', [])
    if len(results) != 1 or len(results[0].get('demand_evaluations', [])) != 1:
        raise ValueError('Unexpected E-Routing result shape')
    evaluation = results[0]['demand_evaluations'][0]
    if request.include_layer1:
        failure = evaluation.get('failure_analysis')
        if not isinstance(failure, dict) or not isinstance(failure.get('tested_links'), int) or failure['tested_links'] <= 0:
            raise ValueError('Layer-1 analysis was requested but no tested-link evidence was produced')
    else:
        failure = None
    links = evaluation.get('link_utilizations')
    if not isinstance(links, list):
        raise ValueError('Updated C++ JSON writer must provide link_utilizations')
    return {
        'scenario_id': request.scenario_id,
        'status': 'completed',
        'execution_mode': 'persisted_routing_scheme',
        'scenario_type': request.scenario_type,
        'demand_model': request.demand_model,
        'demand_scale': request.demand_scale,
        'capacity_upgrades': [{'edge_id': e, 'capacity_factor': f} for e, f in request.capacity_upgrades],
        'evidence': {'result_path': str(result_path), 'routing_table': str(table), 'graph': str(graph)},
        'metrics': {'congestion': evaluation.get('congestion'), 'layer1': failure, 'links': links, 'overloaded_links': [x for x in links if x.get('overloaded')], 'overloaded_count': sum(bool(x.get('overloaded')) for x in links)},
        'limitations': [
            'Fixed persisted routing is evaluated on an intact topology with hypothetical traffic and capacity changes.',
            'Layer-1 reports static exposure, not recovery or rerouting after failure.',
            'Solver provenance is a trusted application label, not authenticated binary metadata.',
            'No live telemetry or network modifications are performed.'
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument('--engine', required=True, type=Path)
    parser.add_argument('--graph', required=True, type=Path)
    parser.add_argument('--routing-table', required=True, type=Path)
    parser.add_argument('--solver', required=True)
    parser.add_argument('--workspace', type=Path, default=Path('results/scenarios'))
    parser.add_argument('--request', type=Path, required=True)
    parser.add_argument('--seed', type=int, default=42)
    args = parser.parse_args()
    request = ScenarioRequest.from_dict(json.loads(args.request.read_text(encoding='utf-8')))
    trusted = TrustedRun(args.engine, args.graph, args.routing_table, args.solver, args.workspace, args.seed)
    print(json.dumps(execute(request, trusted), indent=2, allow_nan=False))


if __name__ == '__main__':
    main()
