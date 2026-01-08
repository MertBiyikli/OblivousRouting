#!/usr/bin/env bash
set -euo pipefail

IMAGE="oblivious-routing:latest"
# NEW: iterate over *all* datasets recursively
DATASET_ROOT="$(pwd)/experiments/datasets"
OUT_DIR="$(pwd)/results/all_datasets"
CSV="$OUT_DIR/results.csv"

SOLVERS_RAW="${1:-electrical,frt,ckr,random_mst,cohen}"
DEMAND="${2:-gravity}"

SOLVERS="$(echo "$SOLVERS_RAW" | tr -d '[:space:]')"

mkdir -p "$OUT_DIR"

# Updated CSV schema (now includes solve/transformation) times separately
# NEW: include dataset + relative graph path for better semantics/debugging
echo "dataset,graph,solver,num_edges,total_time_ms,solve_time_ms,transformation_time_ms,mwu_iterations,avg_oracle_time_ms,achieved_congestion,offline_opt_value" > "$CSV"

# NEW: recursive find (deterministic order)
mapfile -t GRAPHS < <(find "$DATASET_ROOT" -type f -name "*.lgf" | sort)

if [ ${#GRAPHS[@]} -eq 0 ]; then
  echo "ERROR: No .lgf files found under $DATASET_ROOT"
  exit 1
fi

RUN_ID="$(date +%Y%m%d_%H%M%S)"

for g in "${GRAPHS[@]}"; do
  # relative path inside datasets/
  rel_path="${g#$DATASET_ROOT/}"

  # dataset name = first directory component (e.g., Backbone/..., Rocketfuel_Topologies/...)
  dataset="${rel_path%%/*}"

  # filename only (used for log name)
  base="$(basename "$g")"

  # make log filename safe + unique
  safe_rel="${rel_path//\//__}"
  log="$OUT_DIR/${safe_rel%.lgf}_${RUN_ID}.log"

  echo "[RUN] $rel_path"

  # Mount the dataset root once; pass the relative path inside it
  docker run --rm --runtime=runc \
    -e OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}" \
    -v "$DATASET_ROOT":/app/graphs \
    "$IMAGE" \
    "$SOLVERS" "graphs/$rel_path" "$DEMAND" \
    > "$log" 2>&1

  awk -v graph="$base" '
  function flush_row() {
    if (solver != "") {
      printf "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n",
        graph, solver, edges, total_time, solve_time, transformation_time,
        mwu, avg_oracle, achieved, offline
    }
  }

  BEGIN {
    solver="";
    edges="NaN";
    total_time="NaN";
    solve_time="NaN";
    transformation_time="NaN";
    mwu="NaN";
    avg_oracle="NaN";
    achieved="NaN";
    offline="NaN";
  }

  # "Graph loaded: 18 nodes, 66 edges."
  match($0, /^Graph loaded: [0-9]+ nodes, ([0-9]+) edges\./, m) {
    edges = m[1];
    next
  }

  # Start of solver block
  /^=== Running solver:/ {
    flush_row();

    solver = $4;
    total_time="NaN";
    solve_time="NaN";
    transformation_time="NaN";
    mwu="NaN";
    avg_oracle="NaN";
    achieved="NaN";
    offline="NaN";
    next
  }

  match($0, /^Total running time: ([0-9.]+) ms$/, m) { total_time = m[1]; next }
  match($0, /^Solve time: ([0-9.]+) ms$/, m)         { solve_time = m[1]; next }
  match($0, /^Transformation time: ([0-9.]+) ms$/, m){ transformation_time = m[1]; next }
  match($0, /^MWU iterations: ([0-9]+)$/, m)         { mwu = m[1]; next }
  match($0, /^Average oracle time: ([0-9.]+) ms$/, m){ avg_oracle = m[1]; next }

  # "(offline / achieved)"
  match($0, /\(([0-9.]+) \/ ([0-9.]+)\)/, m) {
    offline = m[1];
    achieved  = m[2];
    next
  }

  END { flush_row(); }
  ' "$log" >> "$CSV"

  echo "[DONE] $base"
done

echo "CSV written to: $CSV"
