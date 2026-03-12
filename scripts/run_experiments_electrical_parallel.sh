#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   ./run_experiments_electrical_parallel.sh \
#       --dataset <dir-or-file> \
#       --num_threads "1,2,4,8" \
#       --out results/out.csv

BIN="${BIN:-./build/oblivious_routing}"
SOLVER="electrical_parallel"
DATASET=""
NUM_THREADS_LIST="1"
OUT_CSV="results/results_electrical_parallel.csv"
TIMEOUT="${TIMEOUT:-120m}"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --dataset)               DATASET="${2-}";         shift 2 ;;
    --num_threads|--num-threads) NUM_THREADS_LIST="${2-}"; shift 2 ;;
    --out)                   OUT_CSV="${2-}";          shift 2 ;;
    --bin)                   BIN="${2-}";              shift 2 ;;
    *) echo "Unknown argument: $1"; exit 1 ;;
  esac
done

[[ -z "${DATASET:-}" ]] && { echo "ERROR: --dataset is required."; exit 1; }
[[ ! -x "$BIN" ]]       && { echo "ERROR: binary not found: $BIN"; exit 1; }

# macOS needs gtimeout from coreutils
TBIN="timeout"
command -v "$TBIN" >/dev/null 2>&1 || TBIN="gtimeout"
command -v "$TBIN" >/dev/null 2>&1 || { echo "ERROR: install coreutils (brew install coreutils)"; exit 1; }

# collect .lgf files
GRAPHS=()
if [[ -f "$DATASET" ]]; then
  GRAPHS=("$DATASET")
  DATASET_ROOT="$(cd "$(dirname "$DATASET")" && pwd)"
else
  [[ -d "$DATASET" ]] || { echo "ERROR: dataset not found: $DATASET"; exit 1; }
  DATASET_ROOT="$(cd "$DATASET" && pwd)"
  mapfile -t GRAPHS < <(find "$DATASET_ROOT" -type f -name "*.lgf" | sort)
fi
[[ ${#GRAPHS[@]} -gt 0 ]] || { echo "ERROR: no .lgf files found in $DATASET"; exit 1; }

# parse thread list
IFS=',' read -ra THREADS <<< "$(echo "$NUM_THREADS_LIST" | tr -d ' ')"

mkdir -p "$(dirname "$OUT_CSV")"

echo "dataset,graph,num_nodes,num_edges,num_threads,total_time_ms,solve_time_ms,transformation_time_ms,mwu_iterations,avg_oracle_time_ms,oblivious_ratio,status" > "$OUT_CSV"

for g in "${GRAPHS[@]}"; do
  g_abs="$(cd "$(dirname "$g")" && pwd)/$(basename "$g")"
  graph_name="$(basename "$g_abs" .lgf)"

  for t in "${THREADS[@]}"; do
    echo "[RUN] $graph_name | threads=$t"

    output="$("$TBIN" "$TIMEOUT" "$BIN" "$SOLVER" "$g_abs" "$t" 2>&1)" && status="OK" || {
      ec=$?
      [[ $ec -eq 124 || $ec -eq 137 || $ec -eq 143 ]] && status="TIMEOUT" || status="ERROR_$ec"
    }

    nodes=$(echo "$output"  | awk '/^Graph loaded:/{print $3}')
    edges=$(echo "$output"  | awk '/^Graph loaded:/{print $5}')
    total=$(echo "$output"  | awk '/^Total running time:/{print $4}')
    solve=$(echo "$output"  | awk '/^Solve time:/{print $3}')
    transf=$(echo "$output" | awk '/^Transformation time:/{print $3}')
    mwu=$(echo "$output"    | awk '/^MWU iterations:/{print $3}')
    oracle=$(echo "$output" | awk '/^Average oracle time:/{print $4}')
    ratio=$(echo "$output"  | awk '/^Oblivious ratio of the linear routing scheme:/{print $8}')

    nodes="${nodes:-NaN}"; edges="${edges:-NaN}"; total="${total:-NaN}"
    solve="${solve:-NaN}";  transf="${transf:-NaN}"; mwu="${mwu:-NaN}"
    oracle="${oracle:-NaN}"; ratio="${ratio:-NaN}"

    echo "$graph_name,$graph_name,$nodes,$edges,$t,$total,$solve,$transf,$mwu,$oracle,$ratio,$status" >> "$OUT_CSV"
    echo "[DONE] $graph_name | threads=$t | status=$status"
  done
done

echo "Done. CSV: $OUT_CSV"

