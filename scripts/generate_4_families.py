#!/usr/bin/env python3
"""
Generate 4 additional graph families with ~20 graphs each, edge sizes from 100 to 100,000.
Families: Ring, Clique, Torus2D, Grid3D
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

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
        if result.returncode != 0:
            print(f"  ERROR: {result.stderr[:200]}")
            return False
        print(f"  {result.stdout.strip()}")

        # Add randomized capacities
        if seed is not None:
            random.seed(seed)
        augment_edge_capacities(output_path)
        return True
    except subprocess.TimeoutExpired:
        print(f"  TIMEOUT")
        return False

def augment_edge_capacities(lgf_path):
    """
    Replace uniform edge capacities with random values in range [1, 10].
    """
    lgf_path = Path(lgf_path)
    if not lgf_path.exists():
        return

    try:
        with open(lgf_path, 'r') as f:
            lines = f.readlines()

        in_edges = False
        new_lines = []

        for line in lines:
            if line.strip() == '@edges':
                in_edges = True
                new_lines.append(line)
            elif line.strip().startswith('@'):
                in_edges = False
                new_lines.append(line)
            elif in_edges and line.strip() and not line.startswith('label'):
                parts = line.rstrip('\n').split('\t')
                if len(parts) >= 3:
                    if len(parts) == 4:
                        u, v, label = parts[0], parts[1], parts[2]
                        new_capacity = random.randint(1, 10)
                        new_lines.append(f"{u}\t{v}\t{label}\t{new_capacity}\n")
                    else:
                        new_lines.append(line)
                else:
                    new_lines.append(line)
            else:
                new_lines.append(line)

        with open(lgf_path, 'w') as f:
            f.writelines(new_lines)
    except Exception as e:
        pass

def main():
    output_base = Path('/Users/halilibrahim/Desktop/Thesis/ObliviousRouting/experiments/datasets/synthetic')

    # Target edge counts: ~20 graphs from 100 to 100,000 edges
    # Using geometric progression: roughly 1.5x growth
    target_edges = [100, 150, 225, 338, 506, 759, 1139, 1708, 2562, 3844, 5766, 8649,
                    12973, 19460, 29190, 43785, 65678, 98517, 147776, 221664]

    print("=" * 80)
    print("Generating RING topologies with randomized edge capacities...")
    print("=" * 80)
    # Ring: edges = numToR (circular)
    # To get ~target_edges, set numToR ≈ target_edges
    for i, target in enumerate(target_edges, 1):
        numServerPerToR = max(1, target // 100)
        numToR = max(2, target)
        linkCapacity = 1

        params = {
            "numServerPerToR": numServerPerToR,
            "numToR": numToR,
            "linkCapacity": linkCapacity
        }

        output_path = output_base / 'ring' / f'ring_{i:02d}.lgf'
        print(f"Ring {i:02d}: numToR={numToR}, target≈{target}")
        generate_and_augment_graph('ring', params, output_path, seed=1000 + i)

    print("\n" + "=" * 80)
    print("Generating CLIQUE topologies with randomized edge capacities...")
    print("=" * 80)
    # Clique: edges = numToR * (numToR - 1) / 2
    # To get target_edges: numToR ≈ sqrt(2 * target_edges)
    for i, target in enumerate(target_edges, 1):
        numServerPerToR = max(1, target // 500)
        numToR = max(2, int((2 * target) ** 0.5))
        linkCapacity = 1

        params = {
            "numServerPerToR": numServerPerToR,
            "numToR": numToR,
            "linkCapacity": linkCapacity
        }

        output_path = output_base / 'clique' / f'clique_{i:02d}.lgf'
        print(f"Clique {i:02d}: numToR={numToR}, target≈{target}")
        generate_and_augment_graph('clique', params, output_path, seed=2000 + i)

    print("\n" + "=" * 80)
    print("Generating TORUS2D topologies with randomized edge capacities...")
    print("=" * 80)
    # Torus2D: edges = 2 * numRow * numCol
    # To get target_edges: numRow*numCol ≈ target/2, use numRow=numCol for square
    for i, target in enumerate(target_edges, 1):
        numServerPerToR = max(1, target // 200)
        size = max(2, int((target / 2) ** 0.5))
        numRow = size
        numCol = size
        linkCapacity = 1

        params = {
            "numServerPerToR": numServerPerToR,
            "numRow": numRow,
            "numCol": numCol,
            "linkCapacity": linkCapacity
        }

        output_path = output_base / 'torus2d' / f'torus2d_{i:02d}.lgf'
        print(f"Torus2D {i:02d}: {numRow}x{numCol}, target≈{target}")
        generate_and_augment_graph('torus2d', params, output_path, seed=3000 + i)

    print("\n" + "=" * 80)
    print("Generating GRID3D topologies with randomized edge capacities...")
    print("=" * 80)
    # Grid3D: edges ≈ 3 * numRow * numCol * numLev - (edges on boundaries)
    # Approximate: 3 * n^3 ≈ target, so n ≈ (target/3)^(1/3)
    for i, target in enumerate(target_edges, 1):
        numServerPerToR = max(1, target // 100)
        size = max(2, int((target / 3) ** (1/3)))
        numRow = size
        numCol = size
        numLev = max(2, size // 2)
        linkCapacity = 1

        params = {
            "numServerPerToR": numServerPerToR,
            "numRow": numRow,
            "numCol": numCol,
            "numLev": numLev,
            "linkCapacity": linkCapacity
        }

        output_path = output_base / 'grid3d' / f'grid3d_{i:02d}.lgf'
        print(f"Grid3D {i:02d}: {numRow}x{numCol}x{numLev}, target≈{target}")
        generate_and_augment_graph('grid3d', params, output_path, seed=4000 + i)

    print("\n" + "=" * 80)
    print("Done! Generated 4 new graph families with randomized edge capacities.")
    print("=" * 80)

    # Print summary
    for family in ['ring', 'clique', 'torus2d', 'grid3d']:
        family_dir = output_base / family
        files = sorted(family_dir.glob(f'{family}_*.lgf'))
        print(f"\n{family.upper()} graphs ({len(files)} files):")
        total_size = 0
        for f in files:
            size = f.stat().st_size
            total_size += size
            print(f"  {f.name}: {size/1024:.1f} KB")
        print(f"  Total: {total_size/1024/1024:.1f} MB")

if __name__ == '__main__':
    main()

