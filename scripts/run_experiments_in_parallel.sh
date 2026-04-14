#!/usr/bin/env bash
set -euo pipefail

# Absolute path to the directory containing this script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

DATASET=""
SOLVERS=""
JOBS="${JOBS:-8}"
OUT_CSV="${OUT_CSV:-results/results_parallel.csv}"
TIMEOUT="${TIMEOUT:-120m}"
BIN="${BIN:-./build/oblivious_routing}"

DEMANDS=""
DEMAND_PROVIDED=0

export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-1}"
export NUMEXPR_NUM_THREADS="${NUMEXPR_NUM_THREADS:-1}"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --dataset)      DATASET="${2-}"; shift 2 ;;
    --solvers)      SOLVERS="${2-}"; shift 2 ;;
    --demand)       DEMANDS="${2-}"; DEMAND_PROVIDED=1; shift 2 ;;
    --demands)      DEMANDS="${2-}"; DEMAND_PROVIDED=1; shift 2 ;;
    -j|--jobs)      JOBS="${2-}"; shift 2 ;;
    --out)          OUT_CSV="${2-}"; shift 2 ;;
    --timeout)      TIMEOUT="${2-}"; shift 2 ;;
    --bin)          BIN="${2-}"; shift 2 ;;
    -h|--help)
      echo "Usage:"
      echo "  $0 --dataset <dir-or-file> [--solvers \"frt,ckr\"] [-j 8] [--out results.csv]"
      echo "  $0 --dataset <dir-or-file> [--solvers \"frt,ckr\"] --demands \"gravity,uniform\" [-j 8]"
      echo "Each graph is one job. All solvers and demands are passed in one binary call per graph."
      echo "CSV has one row per (graph × solver × demand_model)."
      exit 0
      ;;
    *) echo "Unknown arg: $1"; exit 1 ;;
  esac
done

if [[ -z "$DATASET" ]]; then
  echo "ERROR: --dataset is required"; exit 1
fi

if [[ ! -x "$BIN" ]]; then
  echo "ERROR: binary not found/executable: $BIN"
  echo "Hint: cmake -S . -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build -j"
  exit 1
fi

# Collect graphs
GRAPHS=()
if [[ -f "$DATASET" ]]; then
  GRAPHS=("$DATASET")
else
  mapfile -t GRAPHS < <(find "$DATASET" -type f -name "*.lgf" | sort)
fi
[[ ${#GRAPHS[@]} -gt 0 ]] || { echo "ERROR: no .lgf found in $DATASET"; exit 1; }

mkdir -p results/parts results/logs

RUN_TS="$(date +%Y%m%d_%H%M%S)"
RUNLIST="results/parts/jobs_${RUN_TS}.txt"
: > "$RUNLIST"

# One job per graph — binary handles all solvers+demands in one call
for g in "${GRAPHS[@]}"; do
  g_abs="$(cd "$(dirname "$g")" && pwd)/$(basename "$g")"
  tag="$(basename "${g_abs%.lgf}")"
  out_part="results/parts/${tag}_${RUN_TS}.csv"
  echo "${g_abs}|${out_part}" >> "$RUNLIST"
done

echo "[INFO] Jobs: $(wc -l < "$RUNLIST")  | parallelism: $JOBS"
echo "[INFO] Joblist: $RUNLIST"

# Build the demand flag to forward (empty string if not provided)
DEMANDS_FLAG=""
if [[ "$DEMAND_PROVIDED" -eq 1 && -n "${DEMANDS:-}" ]]; then
  DEMANDS_FLAG="--demands $(echo "$DEMANDS" | tr -d '[:space:]')"
fi


run_one() {
  local g_abs="$1" out_part="$2"
  # shellcheck disable=SC2086
  BIN="$BIN" TIMEOUT="$TIMEOUT" OUT_DIR="results/logs" \
    "$SCRIPT_DIR/run_experiments_using_binary.sh" \
      --solvers "$SOLVERS" \
      --dataset "$g_abs" \
      --out     "$out_part" \
      --bin     "$BIN" \
      $DEMANDS_FLAG
}
export -f run_one
export BIN TIMEOUT SOLVERS DEMANDS_FLAG SCRIPT_DIR

if command -v parallel >/dev/null 2>&1; then
  parallel -j "$JOBS" --colsep '\|' run_one {1} {2} :::: "$RUNLIST"
else
  # Fallback: xargs
  while IFS='|' read -r g_abs out_part; do
    run_one "$g_abs" "$out_part" &
    # Throttle to $JOBS concurrent jobs
    while [[ $(jobs -r | wc -l) -ge $JOBS ]]; do sleep 0.5; done
  done < "$RUNLIST"
  wait
fi

mkdir -p "$(dirname "$OUT_CSV")"
# Collect only the part files produced by THIS run
THIS_RUN_PARTS=()
while IFS='|' read -r _ out_part; do
  [[ -f "$out_part" ]] && THIS_RUN_PARTS+=("$out_part")
done < "$RUNLIST"

if [[ ${#THIS_RUN_PARTS[@]} -gt 0 ]]; then
  head -n 1 "${THIS_RUN_PARTS[0]}" > "$OUT_CSV"
  for p in "${THIS_RUN_PARTS[@]}"; do
    tail -n +2 "$p" >> "$OUT_CSV"
  done
  echo "[DONE] Merged CSV : $OUT_CSV"
else
  echo "[WARN] No part CSVs found — nothing merged."
fi
echo "[DONE] Logs dir   : results/logs"

