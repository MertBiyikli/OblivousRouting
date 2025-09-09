import os
from pathlib import Path


def parse_and_sort_lgf_files(input_dir, output_dir):
    input_dir = Path(input_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    lgf_files = list(input_dir.rglob("*.lgf"))

    for filepath in lgf_files:
        print(f"Processing {filepath}...")

        try:
            with open(filepath, "r", encoding="utf-8", errors="ignore") as f:
                lines = [line.strip() for line in f if line.strip()]
        except Exception as e:
            print(f"Could not read {filepath}: {e}")
            continue

        vertices = set()
        edges_raw = []

        in_nodes = False
        in_arcs = False
        read_arc_header = False
        arc_header = None
        cap_idx = None

        for line in lines:
            lower = line.lower()
            if not in_nodes and lower.startswith("@nodes"):
                in_nodes = True
                in_arcs = False
                continue
            elif in_nodes and (lower.startswith("@arcs") or lower.startswith("@edges")):
                in_nodes = False
                in_arcs = True
                continue
            elif in_nodes:
                parts = line.split()
                for p in parts[:2]:
                    try:
                        nid = int(p)
                        vertices.add(nid)
                    except:
                        continue

        for line in lines:
            lower = line.lower()
            if not in_arcs and (lower.startswith("@arcs") or lower.startswith("@edges")):
                in_arcs = True
                read_arc_header = False
                continue

            if in_arcs and not read_arc_header:
                arc_header = line
                header_tokens = arc_header.split()
                for idx, tok in enumerate(header_tokens):
                    if "capacity" in tok.lower() or "cost" in tok.lower():
                        cap_idx = idx + 2  # +2: node1 + node2 precede column names
                read_arc_header = True
                continue

            if in_arcs:
                parts = line.split()
                if len(parts) < 2:
                    continue
                try:
                    u = int(parts[0])
                    v = int(parts[1])
                    u, v = min(u, v), max(u, v)
                except:
                    continue

                capacity = 1.0
                if cap_idx is not None and cap_idx < len(parts):
                    try:
                        capacity = float(parts[cap_idx])
                    except:
                        pass

                edges_raw.append((u, v, capacity))

        # Deduplicate edges
        edges_raw.sort(key=lambda x: (x[0], x[1]))
        seen_uv = set()
        unique_edges = []
        for u, v, c in edges_raw:
            if (u, v) not in seen_uv:
                seen_uv.add((u, v))
                unique_edges.append((u, v, c))

        sorted_nodes = sorted((v, v) for v in vertices)
        out_file = output_dir / filepath.name
        write_sorted_lgf(out_file, sorted_nodes, unique_edges, arc_header)
        print(f"  ✓ Wrote sorted file to {out_file}")


def write_sorted_lgf(out_path, sorted_nodes, sorted_edges, arc_header):
    with open(out_path, "w", encoding="utf-8") as f:
        f.write("@nodes\n")
        f.write("label\tid\n")
        for label, nid in sorted_nodes:
            f.write(f"{label}\t{nid}\n")

        f.write("\n@arcs\n")
        if arc_header:
            f.write(arc_header + "\n")
        for u, v, capacity in sorted_edges:
            f.write(f"{u}\t{v}\t{capacity}\n")


# Run the script (adjust path if needed)
parse_and_sort_lgf_files(input_dir=".", output_dir="sorted_lgf")
