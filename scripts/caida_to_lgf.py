#!/usr/bin/env python3
"""
Convert CAIDA AS relationship dataset to LGF (Ligra Graph Format).

CAIDA format (directed):
    # Comment lines starting with #
    FromASN  ToASN  Relationship
    8563     3216   -1
    ...

LGF format (undirected, no duplicates):
    @nodes
    id label
    0  "node_0"
    1  "node_1"
    ...

    @edges
    id  u  v  capacity  distance
    0   0  1  1.0       1.0
    1   2  3  1.0       1.0
    ...

The script converts directed edges to undirected representation.
Each undirected edge (u,v) is stored only once in canonical form (u < v).
The graph parser handles treating these as undirected internally.
"""

import argparse
import sys
from pathlib import Path
from collections import defaultdict


def parse_caida_file(input_file):
    """
    Parse CAIDA AS relationship file.
    Returns: nodes (set), edges (list of tuples)

    Note: Creates only one edge per undirected pair (u < v).
    The graph parser handles treating this as undirected internally.
    """
    nodes = set()
    edges = []
    edge_set = set()  # To avoid duplicate edges

    print(f"Reading CAIDA file: {input_file}")

    with open(input_file, 'r') as f:
        line_count = 0
        for line in f:
            line = line.strip()

            # Skip comments and empty lines
            if not line or line.startswith('#'):
                continue

            line_count += 1
            if line_count % 10000 == 0:
                print(f"  Processed {line_count} edges...")

            # Parse: FromASN ToASN Relationship
            parts = line.split()
            if len(parts) < 2:
                print(f"Warning: Skipping malformed line: {line}")
                continue

            try:
                from_asn = int(parts[0])
                to_asn = int(parts[1])
                # relationship = int(parts[2]) if len(parts) > 2 else 0
            except ValueError:
                print(f"Warning: Skipping line with non-integer AS numbers: {line}")
                continue

            # Add nodes
            nodes.add(from_asn)
            nodes.add(to_asn)

            # Add edge (undirected: store canonical form u < v)
            u, v = (from_asn, to_asn) if from_asn < to_asn else (to_asn, from_asn)
            edge_key = (u, v)

            if edge_key not in edge_set:
                edges.append((u, v))
                edge_set.add(edge_key)

    print(f"Parsed {line_count} directed edges from file")
    print(f"Found {len(nodes)} unique nodes")
    print(f"Created {len(edges)} undirected edges (no duplicates)")

    return nodes, edges


def create_node_mapping(nodes):
    """
    Create mapping from AS number to sequential node ID (0, 1, 2, ...).
    Returns: {as_num -> node_id}
    """
    sorted_nodes = sorted(list(nodes))
    node_map = {asn: idx for idx, asn in enumerate(sorted_nodes)}
    return node_map


def write_lgf_file(output_file, nodes, edges, node_map):
    """
    Write graph in LGF format.
    """
    print(f"\nWriting LGF file: {output_file}")

    with open(output_file, 'w') as f:
        # Write header comment
        f.write("@comment\n")
        f.write(f"# Converted from CAIDA AS relationship dataset\n")
        f.write(f"# Nodes: {len(nodes)}\n")
        f.write(f"# Edges: {len(edges)}\n")
        f.write("@comment\n\n")

        # Write nodes section
        f.write("@nodes\n")
        f.write("id\tlabel\n")
        for asn in sorted(nodes):
            node_id = node_map[asn]
            f.write(f"{node_id}\t\"AS{asn}\"\n")

        f.write("\n")

        # Write edges section
        f.write("@edges\n")
        f.write("id\tu\tv\tcapacity\tdistance\n")

        edge_id = 0
        for u_asn, v_asn in edges:
            u_id = node_map[u_asn]
            v_id = node_map[v_asn]
            capacity = 1.0
            distance = 1.0
            f.write(f"{edge_id}\t{u_id}\t{v_id}\t{capacity}\t{distance}\n")
            edge_id += 1

    print(f"LGF file written successfully!")
    print(f"  Total nodes: {len(nodes)}")
    print(f"  Total edges: {len(edges)}")
    print(f"  Output file: {output_file}")


def main():
    parser = argparse.ArgumentParser(
        description="Convert CAIDA AS relationship dataset to LGF format",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python caida_to_lgf.py -i as-caida20071112.txt -o as-caida20071112.lgf
  python caida_to_lgf.py --input caida.txt --output caida.lgf
        """
    )

    parser.add_argument(
        '-i', '--input',
        required=True,
        type=str,
        help='Input CAIDA AS relationship file'
    )
    parser.add_argument(
        '-o', '--output',
        required=True,
        type=str,
        help='Output LGF file'
    )

    args = parser.parse_args()

    input_file = Path(args.input)
    output_file = Path(args.output)

    # Validate input file exists
    if not input_file.exists():
        print(f"Error: Input file not found: {input_file}")
        sys.exit(1)

    # Create output directory if it doesn't exist
    output_file.parent.mkdir(parents=True, exist_ok=True)

    # Parse CAIDA file
    nodes, edges = parse_caida_file(input_file)

    # Create node mapping
    node_map = create_node_mapping(nodes)

    # Write LGF file
    write_lgf_file(output_file, nodes, edges, node_map)

    print("\n✓ Conversion completed successfully!")


if __name__ == '__main__':
    main()

