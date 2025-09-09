import random
import os
from pathlib import Path

def generate_connected_graph(n, m):
    """Generates a connected undirected graph with n nodes and m edges."""
    assert m >= n - 1, "Need at least n-1 edges for connectivity"
    edges = set()
    # First build a spanning tree
    nodes = list(range(n))
    random.shuffle(nodes)
    for i in range(1, n):
        u = nodes[i]
        v = nodes[random.randint(0, i - 1)]
        if u > v:
            u, v = v, u
        edges.add((u, v))

    # Add random edges
    while len(edges) < m:
        u = random.randint(0, n - 1)
        v = random.randint(0, n - 1)
        if u == v:
            continue
        if u > v:
            u, v = v, u
        edges.add((u, v))

    return list(edges)

def write_lgf_file(filename, n, edges, weights):
    with open(filename, "w") as f:
        f.write("@nodes\n")
        for i in range(n):
            f.write(f"{i}\n")

        f.write("\n@edges\n")
        for (u, v), w in zip(edges, weights):
            f.write(f"{u}\t{v}\t{w}\n")


def main():
    outdir = Path("generated_graphs")
    outdir.mkdir(exist_ok=True)

    sizes = [50, 60, 75, 100, 120, 140, 150, 180, 200]
    for n in sizes:
        m = int(n * 1.5)  # ~1.5 * n edges (moderate density)
        edges = generate_connected_graph(n, m)
        weights = [random.randint(1, 20) for _ in edges]
        out_path = outdir / f"randgraph_{n:05d}.lgf"
        write_lgf_file(out_path, n, edges, weights)
        print(f"✅ Generated {out_path} with {n} nodes, {len(edges)} edges")

if __name__ == "__main__":
    main()
