#!/usr/bin/env python3
import argparse
import random
import networkx as nx
from pathlib import Path

def read_lgf(filename):
    """Read a .lgf file with @nodes and @arcs sections into a directed graph."""
    G = nx.Graph()
    with open(filename, "r") as f:
        lines = f.readlines()

    node_section = False
    arc_section = False
    header_skipped = False

    for line in lines:
        line = line.strip()
        if line.startswith("@nodes"):
            node_section = True
            arc_section = False
            header_skipped = False
            continue
        if line.startswith("@arcs") or line.startswith("@edges"):
            arc_section = True
            node_section = False
            header_skipped = False
            continue
        if line.startswith("@") or not line:
            continue

        if node_section:
            if not header_skipped:  # skip header line
                header_skipped = True
                continue
            parts = line.split()
            if len(parts) >= 1:
                node_id = parts[0]
                G.add_node(node_id)

        elif arc_section:
            if not header_skipped:  # skip header line
                header_skipped = True
                continue
            parts = line.split()
            if len(parts) >= 2:
                u, v = parts[0], parts[1]
                if u < v:  # avoid adding both (u,v) and (v,u) in undirected graph
                    G.add_edge(u, v)



    return G



def duplicate_graph(G, num_copies):
    """Return a new graph with multiple disjoint copies of G."""
    H = nx.Graph()
    for i in range(num_copies):
        mapping = {n: f"{n}_{i}" for n in G.nodes()}
        Gi = nx.relabel_nodes(G, mapping)
        H = nx.disjoint_union(H, Gi)
    return H

def connect_components(H, num_extra_edges=0):
    """Ensure H becomes a single connected component by adding edges between them."""
    # While more than one component exists, connect them
    while True:
        components = list(nx.connected_components(H))
        if len(components) <= 1:
            break  # already connected

        # Pick two different components
        c1, c2 = random.sample(components, 2)
        u = random.choice(list(c1))
        v = random.choice(list(c2))
        H.add_edge(u, v)

    # Add optional extra edges between random nodes from (now connected) graph
    for _ in range(num_extra_edges):
        u, v = random.sample(list(H.nodes()), 2)
        if not H.has_edge(u, v):
            H.add_edge(u, v)

    return H



def write_lgf(H, filename):
    """Write graph back into .lgf format with @nodes and @arcs."""
    with open(filename, "w") as f:
        f.write("@nodes\n")
        f.write("label\tnum\tname\n")
        for n, data in H.nodes(data=True):
            num = data.get("num", n)
            name = data.get("name", n)
            f.write(f"{n}\t{num}\t{name}\n")

        f.write("@arcs\n")
        f.write("\tlabel\tname\tcost\n")
        for u, v, data in H.edges(data=True):
            label = data.get("label", "0")
            name = data.get("name", f"{u}_{v}")
            cost = data.get("cost", 1)
            f.write(f"{u}\t{v}\n")


def process_file(input_file, output_dir, copies, extra_edges):
    """Process a single .lgf file and save the modified version."""
    G = read_lgf(input_file)
    H = duplicate_graph(G, copies)
    H = connect_components(H, extra_edges)

    output_path = Path(output_dir) / (Path(input_file).stem + "_modified.lgf")
    write_lgf(H, output_path)
    print(f"✅ Saved: {output_path}")

def main():
    parser = argparse.ArgumentParser(description="Process multiple .lgf graphs")
    parser.add_argument("inputs", nargs="+", help="Input .lgf files or directories")
    parser.add_argument("--output_dir", required=True, help="Output directory")
    parser.add_argument("--copies", type=int, default=2, help="Number of graph duplicates")
    parser.add_argument("--extra_edges", type=int, default=5, help="Number of random edges to add between components")
    args = parser.parse_args()


    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)


    lgf_files = []
    for inp in args.inputs:
        p = Path(inp)
        if p.is_file() and p.suffix == ".lgf":
            lgf_files.append(p)
        elif p.is_dir():
            lgf_files.extend(p.rglob("*.lgf"))

    if not lgf_files:
        print("⚠️ No .lgf files found.")
        return

    for f in lgf_files:
        process_file(f, output_dir, args.copies, args.extra_edges)


if __name__ == "__main__":
    main()
