#!/usr/bin/env bash
set -euo pipefail

############################################
# Defaults (override via env or flags)
############################################
OUT_CSV="results/results.csv"
OUT_DIR="${OUT_DIR:-results/logs}"
TIMEOUT="${TIMEOUT:-20m}"                 # can be overridden: TIMEOUT=30m ...
IMAGE="${IMAGE:-oblivious-routing:latest}"
SOLVERS=""
DATASET=""
DEMAND="gravity"

# -----------------------------
# Args (robust long-option parser)
# -----------------------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    --solvers)
      SOLVERS="${2-}"
      shift 2
      ;;
    --dataset)
      DATASET="${2-}"
      shift 2
      ;;
    --demand)
      DEMAND="${2-}"
      shift 2
      ;;
    --out)
      OUT_CSV="${2-}"
      shift 2
      ;;
    -h|--help)
      echo "Usage: $0 --solvers \"electrical,ckr\" --dataset <dir-or-file> [--demand gravity] [--out results.csv]"
      exit 0
      ;;
    *)
      echo "Unknown / misplaced argument: $1"
      exit 1
      ;;
  esac
done

# Normalize common mistake: user writes "/experiments/..." meaning repo-relative
if [[ "${DATASET:-}" == /* ]] && [[ "${DATASET:-}" == /experiments/* ]]; then
  DATASET="${DATASET#/}"   # strip leading slash
fi

if [[ -z "${SOLVERS:-}" || -z "${DATASET:-}" ]]; then
  echo "ERROR: --solvers and --dataset are required."
  echo "Try: $0 --solvers \"electrical,ckr\" --dataset experiments/datasets/Backbone --demand gravity --out results/x.csv"
  exit 1
fi


############################################
# Resolve timeout command (macOS: gtimeout)
############################################
TIMEOUT_BIN="${TIMEOUT_BIN:-timeout}"
if ! command -v "$TIMEOUT_BIN" >/dev/null 2>&1; then
  if command -v gtimeout >/dev/null 2>&1; then
    TIMEOUT_BIN="gtimeout"
  else
    echo "ERROR: 'timeout' not found."
    echo "  macOS: brew install coreutils  (then gtimeout exists)"
    echo "  or set TIMEOUT_BIN to a valid timeout binary."
    exit 1
  fi
fi

############################################
# Setup output dirs + CSV header
############################################
mkdir -p "$(dirname "$OUT_CSV")"
mkdir -p "$OUT_DIR"

CSV="$OUT_CSV"

if [[ ! -f "$CSV" ]]; then
  echo "dataset,graph,solver,num_edges,total_time_ms,solve_time_ms,transformation_time_ms,mwu_iterations,avg_oracle_time_ms,achieved_congestion,offline_opt_value,status" > "$CSV"
fi

############################################
# Dataset: accept directory OR single file
############################################
DATASET_PATH="$DATASET"
GRAPHS=()

if [[ -f "$DATASET_PATH" ]]; then
  # user gave a single file
  GRAPHS=("$DATASET_PATH")
  DATASET_ROOT="$(cd "$(dirname "$DATASET_PATH")" && pwd)"
else
  # user gave a directory
  if [[ ! -d "$DATASET_PATH" ]]; then
    echo "ERROR: --dataset path not found: $DATASET_PATH"
    exit 1
  fi
  DATASET_ROOT="$(cd "$DATASET_PATH" && pwd)"
  mapfile -t GRAPHS < <(find "$DATASET_ROOT" -type f -name "*.lgf" | sort)
fi

if [[ ${#GRAPHS[@]} -eq 0 ]]; then
  echo "ERROR: No .lgf files found under $DATASET"
  exit 1
fi

############################################
# Solvers array
############################################
SOLVERS_CSV="$(echo "$SOLVERS" | tr -d '[:space:]')"
IFS=',' read -ra SOLVER_ARR <<< "$SOLVERS_CSV"

RUN_ID="$(date +%Y%m%d_%H%M%S)"

############################################
# Main loop
############################################
for g in "${GRAPHS[@]}"; do
  g_abs="$(cd "$(dirname "$g")" && pwd)/$(basename "$g")"

  # Make a repo-relative-ish path w.r.t DATASET_ROOT (best effort)
  if [[ "$g_abs" == "$DATASET_ROOT/"* ]]; then
    rel_path="${g_abs#$DATASET_ROOT/}"
  else
    rel_path="$(basename "$g_abs")"
  fi

  # Dataset label for CSV: top folder under dataset root, or root folder name if single file
  if [[ "$rel_path" == *"/"* ]]; then
    dataset_label="${rel_path%%/*}"
  else
    dataset_label="$(basename "$DATASET_ROOT")"
  fi

  base="$(basename "$g_abs")"

  for solver in "${SOLVER_ARR[@]}"; do
    solver="$(echo "$solver" | tr -d '[:space:]')"
    if [[ -z "$solver" ]]; then
      continue
    fi

    safe_rel="${rel_path//\//__}"
    log="$OUT_DIR/${safe_rel%.lgf}_${solver}_${RUN_ID}.log"

    echo "[RUN] $rel_path | solver=$solver | demand=$DEMAND | timeout=$TIMEOUT"

    status="OK"
    if "$TIMEOUT_BIN" --signal=SIGTERM --kill-after=30s "$TIMEOUT" \
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
      if [[ "$ec" -eq 124 ]]; then
        status="TIMEOUT"
      else
        status="ERROR_$ec"
      fi
    fi

    # Parse log (BSD-awk compatible: no match(..., m) capture arrays)
    awk -v dataset="$dataset_label" -v graph="$rel_path" -v solver="$solver" -v status="$status" '
    BEGIN {
      edges="NaN"; total_time="NaN"; solve_time="NaN"; transformation_time="NaN";
      mwu="NaN"; avg_oracle="NaN"; achieved="NaN"; offline="NaN";
    }

    # Graph loaded: 18 nodes, 66 edges.
    $0 ~ /^Graph loaded: [0-9]+ nodes, [0-9]+ edges\./ {
      tmp=$0
      sub(/^Graph loaded: [0-9]+ nodes, /, "", tmp)
      sub(/ edges\..*$/, "", tmp)
      edges=tmp
      next
    }

    $0 ~ /^Total running time: [0-9.]+ ms$/ {
      tmp=$0; sub(/^Total running time: /,"",tmp); sub(/ ms$/,"",tmp); total_time=tmp; next
    }
    $0 ~ /^Solve time: [0-9.]+ ms$/ {
      tmp=$0; sub(/^Solve time: /,"",tmp); sub(/ ms$/,"",tmp); solve_time=tmp; next
    }
    $0 ~ /^Transformation time: [0-9.]+ ms$/ {
      tmp=$0; sub(/^Transformation time: /,"",tmp); sub(/ ms$/,"",tmp); transformation_time=tmp; next
    }
    $0 ~ /^MWU iterations: [0-9]+$/ {
      tmp=$0; sub(/^MWU iterations: /,"",tmp); mwu=tmp; next
    }
    $0 ~ /^Average oracle time: [0-9.]+ ms$/ {
      tmp=$0; sub(/^Average oracle time: /,"",tmp); sub(/ ms$/,"",tmp); avg_oracle=tmp; next
    }

    # "(achieved / offline)"
    $0 ~ /\([0-9.]+ \/ [0-9.]+\)/ {
      tmp=$0
      sub(/^.*\(/,"",tmp)
      sub(/\).*$/,"",tmp)
      split(tmp, a, " / ")
      achieved=a[2]
      offline=a[1]
      next
    }

    END {
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
echo "Logs written to: $OUT_DIR"
