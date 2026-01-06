#!/usr/bin/env bash
set -euo pipefail

IMAGE="oblivious-routing:latest"
DATASET_DIR="$(pwd)/experiments/datasets/Rocketfuel_Topologies"
OUT_DIR="$(pwd)/results/Rocketfuel_Topologies"
CSV="$OUT_DIR/results.csv"

SOLVERS_RAW="${1:-electrical, frt, ckr}"
DEMAND="${2:-gravity}"

# normalize solver list (remove whitespace)
SOLVERS="$(echo "$SOLVERS_RAW" | tr -d '[:space:]')"

mkdir -p "$OUT_DIR"

# CSV header expected by plot_experiments.py
echo "solver,num_edges,total_time_ms,mwu_iterations,achieved_congestion,offline_opt_value,avg_oracle_time_ms" > "$CSV"

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

  # num_edges from LGF (common header line: "@edges <m>")
  num_edges=$(grep -E "^@edges" "$g" | awk '{print $2}' | head -n 1)
  num_edges=${num_edges:-NaN}

  # Parse metrics from log
  total_time=$(grep "Total running time:" "$log" | awk '{print $4}')
  avg_oracle=$(grep "Average oracle time:" "$log" | awk '{print $5}')
  mwu=$(grep "MWU iterations:" "$log" | awk '{print $3}')

  total_time=${total_time:-NaN}
  avg_oracle=${avg_oracle:-NaN}
  mwu=${mwu:-NaN}

  ratio_line=$(grep "Ratio off the optimal offline solution" "$log" || true)
  achieved_congestion=$(echo "$ratio_line" | awk -F'[()/]' '{print $2}' | xargs)
  offline_opt_value=$(echo "$ratio_line" | awk -F'[()/]' '{print $3}' | xargs)

  achieved_congestion=${achieved_congestion:-NaN}
  offline_opt_value=${offline_opt_value:-NaN}

  # one row per solver
  for solver in $(echo "$SOLVERS" | tr ',' ' '); do
    echo "$solver,$num_edges,$total_time,$mwu,$achieved_congestion,$offline_opt_value,$avg_oracle" >> "$CSV"
  done

  echo "[DONE] $base"
done

echo "======================================"
echo "All runs finished"
echo "CSV written to: $CSV"
