import argparse
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import re

def extract_threads_from_filename(filename: str) -> int:
    match = re.search(r"threads_(\d+)", filename)
    return int(match.group(1)) if match else -1

def count_vertices_from_lgf(graph_path: Path) -> int:
    """Parses the LGF file and counts number of node entries."""
    if not graph_path.exists():
        return -1
    try:
        with graph_path.open("r") as f:
            for line in f:
                if line.strip().startswith("@nodes"):
                    break
            count = 0
            for line in f:
                line = line.strip()
                if line.startswith("@") or line == "":
                    break
                count += 1
            return count
    except Exception as e:
        print(f"❌ Error reading {graph_path}: {e}")
        return -1

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--csv', nargs='+', required=True, help='List of CSV files')
    parser.add_argument('--outdir', required=True, help='Directory to save the plot')
    parser.add_argument('--solver', default='electrical', help='Filter solver name (default: electrical)')
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    all_data = []

    for file in args.csv:
        path = Path(file)
        if not path.exists():
            print(f"⚠️ File not found: {path}")
            continue

        num_threads = extract_threads_from_filename(path.name)
        df = pd.read_csv(path)
        df = df[df["solver"].str.contains(args.solver, case=False)]
        df["threads"] = num_threads
        df["graph_path"] = df["graph"].apply(lambda p: Path(p))
        df["num_vertices"] = df["graph_path"].apply(count_vertices_from_lgf)
        df["graph_name"] = df["graph_path"].apply(lambda p: p.stem)
        all_data.append(df)

    if not all_data:
        print("❌ No valid data loaded.")
        return

    df_all = pd.concat(all_data)
    df_all = df_all[df_all["num_vertices"] > 0]

    # Compute average wall time per graph per thread count
    grouped = df_all.groupby(["threads", "graph_name", "num_vertices"]).agg({"wall_time_s": "mean"}).reset_index()

    # Create a fixed graph order (by num_vertices, then name)
    graph_order = (
        grouped[["graph_name", "num_vertices"]]
        .drop_duplicates()
        .sort_values(["num_vertices", "graph_name"])
        .reset_index(drop=True)
    )
    graph_order["x_index"] = range(len(graph_order))

    # Merge to assign x-axis index
    grouped = pd.merge(grouped, graph_order, on=["graph_name", "num_vertices"], how="left")

    # Plotting
    plt.figure(figsize=(12, 6))
    for threads, group_df in grouped.groupby("threads"):
        plt.scatter(group_df["x_index"], group_df["wall_time_s"],
                    label=f"{threads} threads", s=20)

    # X-axis labels: graph names (optional: also show |V| if you prefer)
    plt.xticks(graph_order["x_index"], [f"{name}\n({n})" for name, n in zip(graph_order["graph_name"], graph_order["num_vertices"])],
               rotation=45, ha="right", fontsize=8)

    plt.xlabel("Graph instances (name + |V|)")
    plt.ylabel("Avg wall time [s]")
    plt.title("Electrical Flow Runtime vs Graph Instances")
    plt.grid(True, linestyle="--", alpha=0.5)
    plt.legend(title="OMP_NUM_THREADS")
    plt.tight_layout()

    output_file = outdir / "runtime_vs_instance_aligned.png"
    plt.savefig(output_file, dpi=180)
    print(f"✅ Plot saved to: {output_file}")

if __name__ == "__main__":
    main()
