#!/usr/bin/env bash
set -euo pipefail

DEFAULT_SOLVERS="electrical,frt,ckr"
DEFAULT_DEMAND="gravity"

SOLVERS_RAW="${1:-$DEFAULT_SOLVERS}"
DEMAND="${2:-$DEFAULT_DEMAND}"

# Normalize solver list: "electrical, frt" → "electrical,frt"
SOLVERS="$(echo "$SOLVERS_RAW" | tr -d '[:space:]')"

# --------------------------------
# Configuration
# --------------------------------
IMAGE="oblivious-routing"
GRAPH_DIR="$(pwd)/experiments/datasets/Backbone"
RESULT_DIR="$(pwd)/results/raw"
mkdir -p "$RESULT_DIR"

RUN_ID="$(date +%Y%m%d_%H%M%S)"

# --------------------------------
# Sanity checks
# --------------------------------
if [ ! -d "$GRAPH_DIR" ]; then
  echo "ERROR: graphs/ directory not found"
  exit 1
fi

shopt -s nullglob
GRAPHS=("$GRAPH_DIR"/*.lgf)
if [ ${#GRAPHS[@]} -eq 0 ]; then
  echo "ERROR: no .lgf files found in graphs/"
  exit 1
fi
shopt -u nullglob

# --------------------------------
# Run experiments
# --------------------------------
for graph in "${GRAPHS[@]}"; do
  graph_name="$(basename "$graph")"

  echo "----------------------------------------------"
  echo "Graph: $graph_name"
  echo "----------------------------------------------"

  docker run --rm \
    -e OMP_NUM_THREADS=1 \
    -v "$GRAPH_DIR":/app/graphs \
    -v "$RESULT_DIR":/app/results \
    "$IMAGE" \
    "$SOLVERS" "graphs/$graph_name" "$DEMAND" \
    | tee "$RESULT_DIR/run_${RUN_ID}_${graph_name}.log"
done

echo "All experiments finished."
