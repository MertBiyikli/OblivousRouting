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

# CSV header
echo "graph,solvers,demand,total_time_ms,avg_oracle_time_ms,mwu_iterations,oblivious_ratio_percent,ratio_num,ratio_den" > "$CSV"

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

  # -------------------------------
  # Parse metrics from log
  # -------------------------------
  total_time=$(grep "Total running time:" "$log" | awk '{print $4}')
  avg_oracle=$(grep "Average oracle time:" "$log" | awk '{print $5}')
  mwu=$(grep "MWU iterations:" "$log" | awk '{print $3}')

  ratio_line=$(grep "Ratio off the optimal offline solution" "$log" || true)
  ratio_pct=$(echo "$ratio_line" | awk '{print $9}' | tr -d '%')
  ratio_num=$(echo "$ratio_line" | awk -F'[()/]' '{print $2}')
  ratio_den=$(echo "$ratio_line" | awk -F'[()/]' '{print $3}')

  # fallback safety
  total_time=${total_time:-NaN}
  avg_oracle=${avg_oracle:-NaN}
  mwu=${mwu:-NaN}
  ratio_pct=${ratio_pct:-NaN}
  ratio_num=${ratio_num:-NaN}
  ratio_den=${ratio_den:-NaN}

  echo "$base,$SOLVERS,$DEMAND,$total_time,$avg_oracle,$mwu,$ratio_pct,$ratio_num,$ratio_den" >> "$CSV"

  echo "[DONE] $base"
done

echo "======================================"
echo "All runs finished"
echo "CSV written to: $CSV"
