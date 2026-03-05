#!/usr/bin/env bash
set -euo pipefail

DATASET=""
SOLVERS="${SOLVERS:-0,1,2,3,4,6,7}"
DEMANDS="${DEMANDS:-gravity,gaussian,uniform,bimodal}"
JOBS="${JOBS:-8}"
OUT_CSV="${OUT_CSV:-results/results_parallel.csv}"
TIMEOUT="${TIMEOUT:-600m}"
BIN="${BIN:-./build/oblivious_routing}"

# Optional: make dynamic libs visible (OR-Tools, Boost, protobuf in ~/.local)
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-1}"
export NUMEXPR_NUM_THREADS="${NUMEXPR_NUM_THREADS:-1}"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --dataset) DATASET="${2-}"; shift 2 ;;
    --solvers) SOLVERS="${2-}"; shift 2 ;;
    --demands) DEMANDS="${2-}"; shift 2 ;;
    -j|--jobs) JOBS="${2-}"; shift 2 ;;
    --out) OUT_CSV="${2-}"; shift 2 ;;
    --timeout) TIMEOUT="${2-}"; shift 2 ;;
    --bin) BIN="${2-}"; shift 2 ;;
    -h|--help)
      echo "Usage:"
      echo "  $0 --dataset <dir-or-file> [--solvers \"0,1,2\"] [--demands \"gravity,uniform\"] [-j 8] [--out results.csv]"
      exit 0
      ;;
    *) echo "Unknown arg: $1"; exit 1 ;;
  esac
done

if [[ -z "$DATASET" ]]; then
  echo "ERROR: --dataset is required"
  exit 1
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
RUNLIST="results/parts/jobs_$(date +%Y%m%d_%H%M%S).txt"
: > "$RUNLIST"

SOLVERS_CSV="$(echo "$SOLVERS" | tr -d '[:space:]')"
DEMANDS_CSV="$(echo "$DEMANDS" | tr -d '[:space:]')"
IFS=',' read -ra SOLVER_ARR <<< "$SOLVERS_CSV"
IFS=',' read -ra DEMAND_ARR <<< "$DEMANDS_CSV"

# Create job list: solver|graph|demand|outcsv
for g in "${GRAPHS[@]}"; do
  g_abs="$(cd "$(dirname "$g")" && pwd)/$(basename "$g")"
  for s in "${SOLVER_ARR[@]}"; do
    s="$(echo "$s" | tr -d '[:space:]')"; [[ -n "$s" ]] || continue
    for d in "${DEMAND_ARR[@]}"; do
      d="$(echo "$d" | tr -d '[:space:]')"; [[ -n "$d" ]] || continue
      tag="$(basename "${g_abs%.lgf}")__s${s}__${d}"
      out_part="results/parts/${tag}.csv"
      echo "${s}|${g_abs}|${d}|${out_part}" >> "$RUNLIST"
    done
  done
done

echo "[INFO] Jobs: $(wc -l < "$RUNLIST")  | parallelism: $JOBS"
echo "[INFO] Joblist: $RUNLIST"

run_one() {
  local s="$1" g="$2" d="$3" out_part="$4"
  # Reuse your existing script + parsing + logging:
  BIN="$BIN" TIMEOUT="$TIMEOUT" OUT_DIR="results/logs" \
    ./scripts/run_experiments_using_binary.sh --solvers "$s" --dataset "$g" --demand "$d" --out "$out_part"
}
export -f run_one
export BIN TIMEOUT

if command -v parallel >/dev/null 2>&1; then
  parallel -j "$JOBS" --colsep '|' run_one {1} {2} {3} {4} :::: "$RUNLIST"
else
  # fallback without GNU parallel
  cat "$RUNLIST" | xargs -I{} -P "$JOBS" bash -lc '
    IFS="|" read -r s g d out <<< "{}"
    run_one "$s" "$g" "$d" "$out"
  '
fi


mkdir -p "$(dirname "$OUT_CSV")"
first_part="$(ls -1 results/parts/*.csv | head -n 1)"
head -n 1 "$first_part" > "$OUT_CSV"
tail -n +2 -q results/parts/*.csv >> "$OUT_CSV"

echo "[DONE] Merged CSV: $OUT_CSV"
echo "[DONE] Logs dir  : results/logs"