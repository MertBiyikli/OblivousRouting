#!/usr/bin/env bash
set -euo pipefail

IMAGE="oblivious-routing:latest"
DATASET_DIR="$(pwd)/experiments/datasets/Rocketfuel_Topologies"
OUT_DIR="$(pwd)/results/Rocketfuel_Topologies"
CSV="$OUT_DIR/results.csv"

SOLVERS_RAW="${1:-electrical, frt, ckr}"
DEMAND="${2:-gravity}"

# normalize solver list
SOLVERS="$(echo "$SOLVERS_RAW" | tr -d '[:space:]')"

mkdir -p "$OUT_DIR"

# CSV schema expected by plot_experiments.py
echo "graph,solver,num_edges,total_time_ms,mwu_iterations,achieved_congestion,offline_opt_value,avg_oracle_time_ms" > "$CSV"


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

  # ------------------------------------------------------------
  # Parse log: one CSV row per solver block
  # ------------------------------------------------------------
 awk -v graph="$base" '
 BEGIN {
   solver="";
   edges="NaN";
   total_time="NaN";
   avg_oracle="NaN";
   mwu="NaN";
   achieved="NaN";
   offline="NaN";
 }

 # Graph loaded: 18 nodes, 66 edges.
 match($0, /^Graph loaded: [0-9]+ nodes, ([0-9]+) edges\./, m) {
   edges = m[1];
   next
 }

 # Start of solver block
 /^=== Running solver:/ {
   if (solver != "") {
     printf "%s,%s,%s,%s,%s,%s,%s,%s\n",
       graph, solver, edges, total_time, mwu, achieved, offline, avg_oracle
   }
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

 # Total running time: 1671 ms
 match($0, /^Total running time: ([0-9.]+) ms$/, m) {
   total_time = m[1];
   next
 }

  # Solve time: 1671 ms
  match($0, /^Solve time: ([0-9.]+) ms$/, m) {
    solve_time = m[1];
    next
  }

 # MWU iterations: 23
 match($0, /^MWU iterations: ([0-9]+)$/, m) {
   mwu = m[1];
   next
 }

  # Average oracle time: 65.4783 ms
  match($0, /^Average oracle time: ([0-9.]+) ms$/, m) {
    avg_oracle = m[1];
    next
  }


 # Ratio ... (1.84449 / 16.7927)
 match($0, /\(([0-9.]+) \/ ([0-9.]+)\)/, m) {
   offline = m[1];
   achieved = m[2];
   next
 }

 END {
   if (solver != "") {
     printf "%s,%s,%s,%s,%s,%s,%s,%s\n",
       graph, solver, edges, total_time, mwu, achieved, offline, avg_oracle
   }
 }
 ' "$log" >> "$CSV"


  echo "[DONE] $base"
done

echo "======================================"
echo "All runs finished"
echo "CSV written to: $CSV"
