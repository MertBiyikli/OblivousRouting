#!/usr/bin/env python3
"""
Generate larger fat-tree and fat-clique graphs with randomized edge capacities.
"""

import sys
import json
import subprocess
import random
from pathlib import Path

def generate_and_augment_graph(topology_type, params, output_path, seed=None):
    """
    Generate a graph using the topology generator, then add randomized edge weights.
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Generate base graph
    cmd = [
        'python3',
        '/Users/halilibrahim/Desktop/Thesis/ObliviousRouting/scripts/generate_topologies.py',
        '-t', topology_type,
        '-o', str(output_path),
        '-p', json.dumps(params)
    ]

    print(f"Generating {topology_type}: {params}")
    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        print(f"  ERROR: {result.stderr}")
        return False

    print(f"  {result.stdout.strip()}")

    # Add randomized capacities
    if seed is not None:
        random.seed(seed)

    augment_edge_capacities(output_path)
    return True

def augment_edge_capacities(lgf_path):
    """
    Replace uniform edge capacities with random values in range [1, 10].
    """
    lgf_path = Path(lgf_path)

    # Read the file
    with open(lgf_path, 'r') as f:
        lines = f.readlines()

    # Find @edges section and modify capacities
    in_edges = False
    new_lines = []

    for i, line in enumerate(lines):
        if line.strip() == '@edges':
            in_edges = True
            new_lines.append(line)
        elif line.strip().startswith('@'):
            in_edges = False
            new_lines.append(line)
        elif in_edges and line.strip() and not line.startswith('label'):
            # Parse edge line
            parts = line.rstrip('\n').split('\t')
            if len(parts) >= 3:
                # parts[0] = u, parts[1] = v, parts[2] = label (arc id)
                # Last part is capacity
                if len(parts) == 4:
                    # format: u v label cost
                    u, v, label = parts[0], parts[1], parts[2]
                    new_capacity = random.randint(1, 10)
                    new_lines.append(f"{u}\t{v}\t{label}\t{new_capacity}\n")
                else:
                    new_lines.append(line)
            else:
                new_lines.append(line)
        else:
            new_lines.append(line)

    # Write back
    with open(lgf_path, 'w') as f:
        f.writelines(new_lines)

def main():
    output_base = Path('/Users/halilibrahim/Desktop/Thesis/ObliviousRouting/experiments/datasets/synthetic')

    print("=" * 80)
    print("Generating LARGER Fat-Tree topologies with randomized edge capacities...")
    print("=" * 80)

    # Generate larger Fat-Trees
    # Target roughly 10k, 20k, 30k, ..., 100k edges
    # Fat-tree edges ≈ numAgg * switchRadix * switchRadix / 2
    for i in range(1, 11):
        switchRadix = 16 + i * 8  # 24, 32, 40, ..., 96
        numAgg = max(2, i)  # 2, 3, 4, ..., 11

        # Ensure numAgg <= switchRadix
        while numAgg > switchRadix:
            switchRadix += 4

        params = {
            "switchRadix": switchRadix,
            "numAgg": numAgg
        }

        output_path = output_base / 'fattree' / f'fattree_large_{i}.lgf'
        generate_and_augment_graph('fattree_partial', params, output_path, seed=100 + i)

    print("\n" + "=" * 80)
    print("Generating LARGER Fat-Clique topologies with randomized edge capacities...")
    print("=" * 80)

    # Generate larger Fat-Cliques
    # Fat-clique edges = 3 * C(numLocalToR*numSubblock, 2) + similar combinations
    # Use moderate parameters that grow
    for i in range(1, 11):
        numServerPerToR = 2
        numLocalToR = 3 + i
        numSubblock = 3 + i
        numBlock = 2 + i // 2

        params = {
            "numServerPerToR": numServerPerToR,
            "numLocalToR": numLocalToR,
            "numSubblock": numSubblock,
            "numBlock": numBlock
        }

        output_path = output_base / 'fatclique' / f'fatclique_large_{i}.lgf'
        generate_and_augment_graph('fatclique', params, output_path, seed=200 + i)

    print("\n" + "=" * 80)
    print("Done! Generated graphs with randomized edge capacities.")
    print("=" * 80)

    # Print summary
    print("\nFat-Tree files:")
    for f in sorted((output_base / 'fattree').glob('fattree_large_*.lgf')):
        size = f.stat().st_size / 1024  # KB
        print(f"  {f.name}: {size:.1f} KB")

    print("\nFat-Clique files:")
    for f in sorted((output_base / 'fatclique').glob('fatclique_large_*.lgf')):
        size = f.stat().st_size / 1024  # KB
        print(f"  {f.name}: {size:.1f} KB")

if __name__ == '__main__':
    main()

