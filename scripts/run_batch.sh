#!/usr/bin/env bash
set -euo pipefail

IMAGE="oblivious-routing:latest"
DATASET_DIR="$(pwd)/experiments/datasets/Backbone"
OUT_DIR="$(pwd)/results/backbone"
SOLVERS_RAW="${1:-electrical, frt, ckr}"
DEMAND="${2:-gravity}"

# normalize solver list: remove whitespace -> "electrical,frt,ckr"
SOLVERS="$(echo "$SOLVERS_RAW" | tr -d '[:space:]')"

mkdir -p "$OUT_DIR"

# collect graphs
shopt -s nullglob
GRAPHS=("$DATASET_DIR"/*.lgf)
shopt -u nullglob

if [ ${#GRAPHS[@]} -eq 0 ]; then
  echo "ERROR: No .lgf files found in $DATASET_DIR"
  exit 1
fi

RUN_ID="$(date +%Y%m%d_%H%M%S)"
MASTER_LOG="$OUT_DIR/run_${RUN_ID}_MASTER.log"

echo "Dataset: $DATASET_DIR" | tee "$MASTER_LOG"
echo "Solvers: $SOLVERS"     | tee -a "$MASTER_LOG"
echo "Demand:  $DEMAND"      | tee -a "$MASTER_LOG"
echo "Graphs:  ${#GRAPHS[@]}"| tee -a "$MASTER_LOG"
echo "----------------------------------------" | tee -a "$MASTER_LOG"

for g in "${GRAPHS[@]}"; do
  base="$(basename "$g")"
  log="$OUT_DIR/${base%.lgf}__${SOLVERS}__${DEMAND}__${RUN_ID}.log"

  echo "[RUN] $base" | tee -a "$MASTER_LOG"

  docker run --rm --runtime=runc \
    -e OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}" \
    -v "$DATASET_DIR":/app/graphs \
    "$IMAGE" \
    "$SOLVERS" "graphs/$base" "$DEMAND" \
    > "$log" 2>&1

  echo "[DONE] $base -> $log" | tee -a "$MASTER_LOG"
done

echo "All runs finished. Master log: $MASTER_LOG"
