#!/usr/bin/env python3
"""
Convert CAIDA .txt files to LGF format with undirected edges.
- Converts all directed edges to undirected edges
- If two directed edges have the same endpoints (u,v) and (v,u), only one undirected edge is created
- Outputs in LEMON LGF format
"""

import os
import sys
from pathlib import Path
from collections import defaultdict

def read_caida_txt(filepath):
    """
    Read CAIDA .txt file and extract edges.
    Returns:
    - node_mapping: dict mapping original node IDs to sequential IDs starting from 0
    - edges: set of edges (as tuples with u <= v for undirected representation)
    """
    edges = set()
    nodes = set()

    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            # Skip comments and empty lines
            if not line or line.startswith('#'):
                continue

            # Parse edge: FromNodeId ToNodeId Relationship
            parts = line.split()
            if len(parts) >= 2:
                try:
                    u = int(parts[0])
                    v = int(parts[1])

                    nodes.add(u)
                    nodes.add(v)

                    # Store undirected edge (canonical form: smaller node first)
                    edge = tuple(sorted([u, v]))
                    edges.add(edge)
                except ValueError:
                    continue

    # Create mapping from original node IDs to sequential IDs starting from 0
    sorted_nodes = sorted(nodes)
    node_mapping = {original_id: new_id for new_id, original_id in enumerate(sorted_nodes)}

    # Remap edges using the new node IDs
    remapped_edges = set()
    for u, v in edges:
        new_u = node_mapping[u]
        new_v = node_mapping[v]
        remapped_edges.add(tuple(sorted([new_u, new_v])))

    return node_mapping, remapped_edges

def write_lgf(filepath, node_count, edges):
    """
    Write nodes and edges to LGF format.
    Each undirected edge is represented as two directed arcs in LGF.
    """
    with open(filepath, 'w') as f:
        # Write nodes section
        f.write("@nodes\n")
        f.write("label\n")
        for node_id in range(node_count):
            f.write(f"{node_id}\n")

        # Write arcs section (each undirected edge becomes two directed arcs)
        f.write("@arcs\n")
        f.write("\tlabel\n")

        arc_id = 0
        for u, v in sorted(edges):
            # Arc u -> v
            f.write(f"{u}\t{v}\t{arc_id}\n")
            arc_id += 1
            # Arc v -> u
            f.write(f"{v}\t{u}\t{arc_id}\n")
            arc_id += 1

def convert_file(input_file, output_file):
    """Convert a single CAIDA .txt file to LGF format"""
    print(f"Converting {input_file}...")

    try:
        node_mapping, edges = read_caida_txt(input_file)
        node_count = len(node_mapping)
        print(f"  Found {node_count} nodes and {len(edges)} undirected edges")

        write_lgf(output_file, node_count, edges)
        print(f"  Wrote to {output_file}")
        return True
    except Exception as e:
        print(f"  Error: {e}")
        return False

def main():
    # Find all .txt files in the caida/as-caida directory
    caida_dir = Path("/Users/halilibrahim/Desktop/Thesis/ObliviousRouting/experiments/datasets/caida/as-caida")
    output_dir = Path("/Users/halilibrahim/Desktop/Thesis/ObliviousRouting/experiments/datasets/caida_lgf/as-caida")

    # Create output directory if it doesn't exist
    output_dir.mkdir(parents=True, exist_ok=True)

    # Find all .txt files
    txt_files = sorted(caida_dir.glob("*.txt"))

    if not txt_files:
        print(f"No .txt files found in {caida_dir}")
        return

    print(f"Found {len(txt_files)} .txt files to convert\n")

    success_count = 0
    for txt_file in txt_files:
        # Generate output filename
        output_file = output_dir / txt_file.stem.replace("as-caida", "as-caida-undirected")
        output_file = output_file.with_suffix('.lgf')

        if convert_file(txt_file, output_file):
            success_count += 1

    print(f"\nConversion complete: {success_count}/{len(txt_files)} files converted successfully")

if __name__ == "__main__":
    main()

