#!/usr/bin/env bash
set -euo pipefail

IMAGE="oblivious-routing:latest"
DATASET_DIR="$(pwd)/experiments/datasets/Rocketfuel_Topologies"   # <-- adjust if needed
OUT_DIR="$(pwd)/results/Rocketfuel_Topologies"
CSV="$OUT_DIR/results.csv"

SOLVERS_RAW="${1:-electrical,frt,ckr,cohen}"
DEMAND="${2:-gravity}"

SOLVERS="$(echo "$SOLVERS_RAW" | tr -d '[:space:]')"

mkdir -p "$OUT_DIR"

# Updated CSV schema (now includes solve/transformation)
echo "graph,solver,num_edges,total_time_ms,solve_time_ms,transformation_time_ms,mwu_iterations,avg_oracle_time_ms,achieved_congestion,offline_opt_value" > "$CSV"

shopt -s nullglob
GRAPHS=("$DATASET_DIR"/*.lgf)
shopt -u nullglob

if [ ${#GRAPHS[@]} -eq 0 ]; then
  echo "ERROR: No .lgf files found in $DATASET_DIR"
  exit 1
fi

RUN_ID="$(date +%Y%m%d_%H%M%S)"

for g in "${GRAPHS[@]}"; do
  base="$(basename "$g")"
  log="$OUT_DIR/${base%.lgf}_${RUN_ID}.log"

  echo "[RUN] $base"

  docker run --rm --runtime=runc \
    -e OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}" \
    -v "$DATASET_DIR":/app/graphs \
    "$IMAGE" \
    "$SOLVERS" "graphs/$base" "$DEMAND" \
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
