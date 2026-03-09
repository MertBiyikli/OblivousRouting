import os
import json
import argparse
import statistics
import networkx as nx
from pathlib import Path


def load_lgf(filepath: str) -> nx.Graph:
    """Parse a LEMON Graph Format (.lgf) file into a NetworkX graph."""
    g = nx.Graph()
    section = None

    with open(filepath, "r") as f:
        for line in f:
            line = line.strip()

            if not line or line.startswith("#"):
                continue

            if line.startswith("@nodes"):
                section = "nodes"
                header_skipped = False
                continue
            elif line.startswith("@edges") or line.startswith("@arcs"):
                section = "edges"
                header_skipped = False
                continue
            elif line.startswith("@"):
                section = None
                continue

            if section == "nodes":
                if not header_skipped:
                    header_skipped = True
                    continue
                parts = line.split()
                if parts:
                    g.add_node(parts[0])

            elif section == "edges":
                if not header_skipped:
                    header_skipped = True
                    continue
                parts = line.split()
                if len(parts) >= 2:
                    g.add_edge(parts[0], parts[1])

    return g


def load_graph(filepath: str) -> nx.Graph:
    """Load a graph from various file formats."""
    ext = Path(filepath).suffix.lower()

    if ext == ".lgf":
        return load_lgf(filepath)
    elif ext == ".edgelist" or ext == ".txt":
        return nx.read_edgelist(filepath)
    elif ext == ".gml":
        return nx.read_gml(filepath)
    elif ext == ".graphml":
        return nx.read_graphml(filepath)
    elif ext == ".gexf":
        return nx.read_gexf(filepath)
    elif ext == ".json":
        with open(filepath, "r") as f:
            data = json.load(f)
        return nx.node_link_graph(data)
    elif ext == ".csv":
        return nx.read_edgelist(filepath, delimiter=",")
    else:
        raise ValueError(f"Unsupported file format: {ext}")


def degree_stats(degrees: list[int]) -> dict:
    """Compute degree statistics from a list of degree values."""
    return {
        "min": min(degrees),
        "max": max(degrees),
        "average": round(statistics.mean(degrees), 2),
        "median": round(statistics.median(degrees), 2),
        "stdev": round(statistics.stdev(degrees), 2) if len(degrees) > 1 else 0.0,
    }


def compute_metrics(graphs: dict[str, nx.Graph]) -> dict:
    """Compute min/max/average number of nodes, edges and degree across all graphs."""
    node_counts = [g.number_of_nodes() for g in graphs.values()]
    edge_counts = [g.number_of_edges() for g in graphs.values()]

    # Per-graph average degree, and all degrees pooled across all graphs
    per_graph_avg_degrees = []
    all_degrees = []
    for g in graphs.values():
        degrees = [d for _, d in g.degree()]
        if degrees:
            per_graph_avg_degrees.append(statistics.mean(degrees))
            all_degrees.extend(degrees)

    metrics = {
        "num_graphs": len(graphs),
        "nodes": {
            "min": min(node_counts),
            "max": max(node_counts),
            "average": round(statistics.mean(node_counts), 2),
            "median": round(statistics.median(node_counts), 2),
            "stdev": round(statistics.stdev(node_counts), 2) if len(node_counts) > 1 else 0.0,
        },
        "edges": {
            "min": min(edge_counts),
            "max": max(edge_counts),
            "average": round(statistics.mean(edge_counts), 2),
            "median": round(statistics.median(edge_counts), 2),
            "stdev": round(statistics.stdev(edge_counts), 2) if len(edge_counts) > 1 else 0.0,
        },
        "degree_pooled": degree_stats(all_degrees) if all_degrees else {},
        "degree_per_graph_avg": {
            "min": round(min(per_graph_avg_degrees), 2),
            "max": round(max(per_graph_avg_degrees), 2),
            "average": round(statistics.mean(per_graph_avg_degrees), 2),
            "median": round(statistics.median(per_graph_avg_degrees), 2),
            "stdev": round(statistics.stdev(per_graph_avg_degrees), 2) if len(per_graph_avg_degrees) > 1 else 0.0,
        } if per_graph_avg_degrees else {},
    }
    return metrics


def print_metrics(dataset_name: str, metrics: dict, graphs: dict[str, nx.Graph]):
    """Pretty-print the computed metrics."""
    print(f"\n{'=' * 60}")
    print(f"  Dataset: {dataset_name}")
    print(f"  Total graphs: {metrics['num_graphs']}")
    print(f"{'=' * 60}")

    print(f"\n{'--- Nodes ' + '-' * 50}")
    print(f"  Min     : {metrics['nodes']['min']}")
    print(f"  Max     : {metrics['nodes']['max']}")
    print(f"  Average : {metrics['nodes']['average']}")
    print(f"  Median  : {metrics['nodes']['median']}")
    print(f"  Std Dev : {metrics['nodes']['stdev']}")

    print(f"\n{'--- Edges ' + '-' * 50}")
    print(f"  Min     : {metrics['edges']['min']}")
    print(f"  Max     : {metrics['edges']['max']}")
    print(f"  Average : {metrics['edges']['average']}")
    print(f"  Median  : {metrics['edges']['median']}")
    print(f"  Std Dev : {metrics['edges']['stdev']}")

    if metrics["degree_pooled"]:
        print(f"\n{'--- Degree (pooled across all nodes & graphs) ' + '-' * 14}")
        print(f"  Min     : {metrics['degree_pooled']['min']}")
        print(f"  Max     : {metrics['degree_pooled']['max']}")
        print(f"  Average : {metrics['degree_pooled']['average']}")
        print(f"  Median  : {metrics['degree_pooled']['median']}")
        print(f"  Std Dev : {metrics['degree_pooled']['stdev']}")

    if metrics["degree_per_graph_avg"]:
        print(f"\n{'--- Degree (per-graph averages) ' + '-' * 28}")
        print(f"  Min     : {metrics['degree_per_graph_avg']['min']}")
        print(f"  Max     : {metrics['degree_per_graph_avg']['max']}")
        print(f"  Average : {metrics['degree_per_graph_avg']['average']}")
        print(f"  Median  : {metrics['degree_per_graph_avg']['median']}")
        print(f"  Std Dev : {metrics['degree_per_graph_avg']['stdev']}")



def analyze_dataset(dataset_path: str):
    """Iterate over all graph files in a dataset directory and compute metrics."""
    dataset_path = Path(dataset_path)

    if not dataset_path.exists():
        raise FileNotFoundError(f"Dataset path not found: {dataset_path}")

    supported_extensions = {".lgf", ".edgelist", ".txt", ".gml", ".graphml", ".gexf", ".json", ".csv"}
    graph_files = [
        f for f in sorted(dataset_path.rglob("*"))
        if f.is_file() and f.suffix.lower() in supported_extensions
    ]

    if not graph_files:
        print(f"No supported graph files found in: {dataset_path}")
        return

    graphs = {}
    skipped = []

    for filepath in graph_files:
        try:
            g = load_graph(str(filepath))
            graphs[filepath.name] = g
        except Exception as e:
            skipped.append(filepath.name)
            print(f"  [SKIP] {filepath.name}: {e}")


    metrics = compute_metrics(graphs)
    print_metrics(dataset_path.name, metrics, graphs)


def main():
    parser = argparse.ArgumentParser(
        description="Analyze graph datasets and report node/edge/degree statistics."
    )
    parser.add_argument(
        "dataset",
        nargs="+",
        help="Path(s) to dataset directory or directories containing graph files.",
    )
    args = parser.parse_args()

    for dataset_path in args.dataset:
        print(f"\nAnalyzing: {dataset_path}")
        analyze_dataset(dataset_path)


if __name__ == "__main__":
    main()
