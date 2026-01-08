#!/usr/bin/env bash
set -euo pipefail

IMAGE="oblivious-routing:latest"

# iterate over *all* datasets recursively
DATASET_ROOT="$(pwd)/experiments/datasets"
OUT_DIR="$(pwd)/results/all_datasets"
CSV="$OUT_DIR/results.csv"

# Default solvers and demand
SOLVERS_RAW="${1:-electrical,frt,ckr,random_mst,cohen}"
DEMAND="${2:-gravity}"

# Timeout per (solver, graph). Override via: TIMEOUT=10m ./collect_data.sh ...
TIMEOUT="${TIMEOUT:-20m}"

# Normalize solver list (keep commas, remove spaces)
SOLVERS_CSV="$(echo "$SOLVERS_RAW" | tr -d '[:space:]')"

mkdir -p "$OUT_DIR"

# CSV schema (adds status)
echo "dataset,graph,solver,num_edges,total_time_ms,solve_time_ms,transformation_time_ms,mwu_iterations,avg_oracle_time_ms,achieved_congestion,offline_opt_value,status" > "$CSV"

# recursive find (deterministic order)
mapfile -t GRAPHS < <(find "$DATASET_ROOT" -type f -name "*.lgf" | sort)
if [ ${#GRAPHS[@]} -eq 0 ]; then
  echo "ERROR: No .lgf files found under $DATASET_ROOT"
  exit 1
fi

RUN_ID="$(date +%Y%m%d_%H%M%S)"

# Split solvers into array
IFS=',' read -ra SOLVER_ARR <<< "$SOLVERS_CSV"

for g in "${GRAPHS[@]}"; do
  rel_path="${g#$DATASET_ROOT/}"
  dataset="${rel_path%%/*}"
  base="$(basename "$g")"

  for solver in "${SOLVER_ARR[@]}"; do
    solver="$(echo "$solver" | tr -d '[:space:]')"
    safe_rel="${rel_path//\//__}"
    log="$OUT_DIR/${safe_rel%.lgf}_${solver}_${RUN_ID}.log"

    echo "[RUN] $rel_path | solver=$solver | timeout=$TIMEOUT"

    # Run single solver with timeout. If timeout triggers, exit code is 124.
    status="OK"
    if timeout --signal=SIGTERM --kill-after=30s "$TIMEOUT" \
      docker run --rm --runtime=runc \
        -e OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}" \
        -v "$DATASET_ROOT":/app/graphs \
        "$IMAGE" \
        "$solver" "graphs/$rel_path" "$DEMAND" \
      > "$log" 2>&1
    then
      status="OK"
    else
      ec=$?
      if [ "$ec" -eq 124 ]; then
        status="TIMEOUT"
      else
        status="ERROR_$ec"
      fi
    fi

    # Parse the log (one solver per log -> simpler & robust)
    awk -v dataset="$dataset" -v graph="$rel_path" -v solver="$solver" -v status="$status" '
    BEGIN {
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
    match($0, /^Graph loaded: [0-9]+ nodes, ([0-9]+) edges\./, m) { edges=m[1]; next }

    match($0, /^Total running time: ([0-9.]+) ms$/, m)        { total_time=m[1]; next }
    match($0, /^Solve time: ([0-9.]+) ms$/, m)                { solve_time=m[1]; next }
    match($0, /^Transformation time: ([0-9.]+) ms$/, m)       { transformation_time=m[1]; next }
    match($0, /^MWU iterations: ([0-9]+)$/, m)                { mwu=m[1]; next }
    match($0, /^Average oracle time: ([0-9.]+) ms$/, m)       { avg_oracle=m[1]; next }

    # IMPORTANT: output is "(achieved / offline)"
    match($0, /\(([0-9.]+) \/ ([0-9.]+)\)/, m) { achieved=m[2]; offline=m[1]; next }

    END {
      # If TIMEOUT/ERROR, keep NaNs for metrics but still write row
      if (status != "OK") {
        total_time="NaN"; solve_time="NaN"; transformation_time="NaN";
        mwu="NaN"; avg_oracle="NaN"; achieved="NaN"; offline="NaN";
      }
      printf "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n",
        dataset, graph, solver, edges, total_time, solve_time, transformation_time,
        mwu, avg_oracle, achieved, offline, status
    }
    ' "$log" >> "$CSV"

    echo "[DONE] $base | $solver | $status"
  done
done

echo "CSV written to: $CSV"
