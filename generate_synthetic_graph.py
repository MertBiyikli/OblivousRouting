#!/usr/bin/env python3
"""
generate_datacenter_topologies.py

Generates synthetic data-center graphs (Fat-Tree, Folded-Clos) and writes .lgf files.

Example:
  python generate_datacenter_topologies.py fat_tree --k 48 --out experiments/datasets/synth_dc
  python generate_datacenter_topologies.py clos --leaves 256 --spines 64 --hosts-per-leaf 48 --out experiments/datasets/synth_dc

Notes:
- Output is an undirected LGF using @nodes and @edges.
- Each edge has a "capacity" attribute.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import List, Tuple, Dict


@dataclass(frozen=True)
class Edge:
    u: int
    v: int
    cap: float


def write_lgf(path: Path, n: int, edges: List[Edge]) -> None:
    """
    Write an undirected LGF file with nodes labeled 0..n-1 and edges with capacities.
    Format is compatible with many LEMON-style readers:
      @nodes
      label
      0
      1
      ...
      @edges
      u v label capacity
      0 1 e1 1
      ...
    """
    path.parent.mkdir(parents=True, exist_ok=True)

    # de-duplicate (u,v) ignoring orientation
    seen = set()
    deduped: List[Edge] = []
    for e in edges:
        a, b = (e.u, e.v) if e.u <= e.v else (e.v, e.u)
        key = (a, b)
        if a == b:
            continue
        if key in seen:
            continue
        seen.add(key)
        deduped.append(Edge(a, b, e.cap))

    with path.open("w", encoding="utf-8") as f:
        f.write("@nodes\n")
        f.write("label\n")
        for i in range(n):
            f.write(f"{i}\n")

        f.write("\n@edges\n")
        f.write("label u v capacity\n")
        for i, e in enumerate(deduped):
            f.write(f"{e.u} {e.v} {i} {e.cap}\n")

    print(f"[WROTE] {path} | nodes={n:,}, edges={len(deduped):,}")


def make_fat_tree(k: int, cap_core_agg: float, cap_agg_edge: float, cap_edge_host: float) -> Tuple[int, List[Edge], Dict[str, Tuple[int, int]]]:
    """
    Standard k-ary Fat-Tree:
      - k pods
      - Each pod: k/2 aggregation + k/2 edge switches
      - Core: (k/2)^2 switches
      - Each edge switch connects to k/2 hosts
    Total hosts: k^3 / 4

    Returns:
      (n, edges, ranges) where ranges maps component name -> [start,end) node id ranges.
    """
    if k <= 0 or k % 2 != 0:
        raise ValueError("Fat-tree requires even k (e.g., 4, 8, 16, 48).")

    p = k
    half = k // 2

    num_core = half * half
    num_agg = p * half
    num_edge = p * half
    num_host = p * half * half

    # Node id layout
    core_start = 0
    agg_start = core_start + num_core
    edge_start = agg_start + num_agg
    host_start = edge_start + num_edge
    n = host_start + num_host

    ranges = {
        "core": (core_start, agg_start),
        "agg": (agg_start, edge_start),
        "edge": (edge_start, host_start),
        "host": (host_start, n),
    }

    def core_id(group: int, idx: int) -> int:
        # core groups: half groups, each of size half
        # group in [0..half-1], idx in [0..half-1]
        return core_start + group * half + idx

    def agg_id(pod: int, a: int) -> int:
        # aggregation switches in each pod
        return agg_start + pod * half + a

    def edge_id(pod: int, e: int) -> int:
        return edge_start + pod * half + e

    def host_id(pod: int, edge_sw: int, h: int) -> int:
        # hosts per edge switch: half
        # host index within pod: edge_sw*half + h
        return host_start + pod * (half * half) + edge_sw * half + h

    edges: List[Edge] = []

    # Core <-> Agg:
    # Each agg switch (pod, a) connects to exactly one switch in each core group.
    # Standard fat-tree wiring: agg index a chooses the core switch within each group.
    # For each pod, each agg a connects to core switches core(group=g, idx=a) for all g.
    for pod in range(p):
        for a in range(half):
            u = agg_id(pod, a)
            for g in range(half):
                v = core_id(g, a)
                edges.append(Edge(u, v, cap_core_agg))

    # Agg <-> Edge within each pod (complete bipartite K_{half,half})
    for pod in range(p):
        for a in range(half):
            u = agg_id(pod, a)
            for e in range(half):
                v = edge_id(pod, e)
                edges.append(Edge(u, v, cap_agg_edge))

    # Edge <-> Hosts
    for pod in range(p):
        for e in range(half):
            u = edge_id(pod, e)
            for h in range(half):
                v = host_id(pod, e, h)
                edges.append(Edge(u, v, cap_edge_host))

    return n, edges, ranges


def make_folded_clos(leaves: int, spines: int, hosts_per_leaf: int,
                     cap_leaf_spine: float, cap_leaf_host: float) -> Tuple[int, List[Edge], Dict[str, Tuple[int, int]]]:
    """
    Two-tier folded-Clos a.k.a. leaf-spine:
      - leaves L
      - spines S
      - each leaf connected to all spines (full bipartite)
      - each leaf connected to H hosts

    Returns (n, edges, ranges).
    """
    if leaves <= 0 or spines <= 0 or hosts_per_leaf < 0:
        raise ValueError("leaves, spines must be >0 and hosts_per_leaf >=0")

    spine_start = 0
    leaf_start = spine_start + spines
    host_start = leaf_start + leaves
    num_hosts = leaves * hosts_per_leaf
    n = host_start + num_hosts

    ranges = {
        "spine": (spine_start, leaf_start),
        "leaf": (leaf_start, host_start),
        "host": (host_start, n),
    }

    def spine_id(s: int) -> int:
        return spine_start + s

    def leaf_id(l: int) -> int:
        return leaf_start + l

    def host_id(l: int, h: int) -> int:
        return host_start + l * hosts_per_leaf + h

    edges: List[Edge] = []

    # Leaf <-> Spine (full bipartite)
    for l in range(leaves):
        u = leaf_id(l)
        for s in range(spines):
            v = spine_id(s)
            edges.append(Edge(u, v, cap_leaf_spine))

    # Leaf <-> Hosts
    for l in range(leaves):
        u = leaf_id(l)
        for h in range(hosts_per_leaf):
            v = host_id(l, h)
            edges.append(Edge(u, v, cap_leaf_host))

    return n, edges, ranges


def main() -> None:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)

    # Fat-tree
    ap_ft = sub.add_parser("fat_tree", help="Generate k-ary fat-tree")
    ap_ft.add_argument("--k", type=int, required=True, help="even parameter k (e.g., 8,16,48)")
    ap_ft.add_argument("--cap-core-agg", type=float, default=1.0)
    ap_ft.add_argument("--cap-agg-edge", type=float, default=1.0)
    ap_ft.add_argument("--cap-edge-host", type=float, default=1.0)
    ap_ft.add_argument("--out", type=Path, required=True, help="output directory")
    ap_ft.add_argument("--name", type=str, default=None, help="filename (default: fat_tree_k{K}.lgf)")

    # Folded Clos (leaf-spine)
    ap_cl = sub.add_parser("clos", help="Generate folded-Clos / leaf-spine")
    ap_cl.add_argument("--leaves", type=int, required=True)
    ap_cl.add_argument("--spines", type=int, required=True)
    ap_cl.add_argument("--hosts-per-leaf", type=int, required=True)
    ap_cl.add_argument("--cap-leaf-spine", type=float, default=1.0)
    ap_cl.add_argument("--cap-leaf-host", type=float, default=1.0)
    ap_cl.add_argument("--out", type=Path, required=True, help="output directory")
    ap_cl.add_argument("--name", type=str, default=None, help="filename (default: clos_L{L}_S{S}_H{H}.lgf)")

    args = ap.parse_args()

    if args.cmd == "fat_tree":
        n, edges, ranges = make_fat_tree(
            k=args.k,
            cap_core_agg=args.cap_core_agg,
            cap_agg_edge=args.cap_agg_edge,
            cap_edge_host=args.cap_edge_host,
        )
        fname = args.name or f"fat_tree_k{args.k}.lgf"
        out_path = args.out / fname
        write_lgf(out_path, n, edges)
        print("[RANGES]", ranges)

    elif args.cmd == "clos":
        n, edges, ranges = make_folded_clos(
            leaves=args.leaves,
            spines=args.spines,
            hosts_per_leaf=args.hosts_per_leaf,
            cap_leaf_spine=args.cap_leaf_spine,
            cap_leaf_host=args.cap_leaf_host,
        )
        fname = args.name or f"clos_L{args.leaves}_S{args.spines}_H{args.hosts_per_leaf}.lgf"
        out_path = args.out / fname
        write_lgf(out_path, n, edges)
        print("[RANGES]", ranges)


if __name__ == "__main__":
    main()
