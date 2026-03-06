#!/usr/bin/env bash
set -euo pipefail

OUT_CSV="results/results.csv"
OUT_DIR="${OUT_DIR:-results/logs}"
TIMEOUT="${TIMEOUT:-120m}"
IMAGE="${IMAGE:-oblivious-routing:latest}"

SOLVERS=""
DATASET=""
DEMANDS=""
DEMAND_PROVIDED=0

ALL_SOLVERS="${ALL_SOLVERS:-frt,ckr,mst,frt_mendel,ckr_mendel}"
ALL_DEMANDS="${ALL_DEMANDS:-gravity,gaussian,uniform,bimodal}"

RUN_ALL=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --all|-all)      RUN_ALL=1; shift ;;
    --solvers)       SOLVERS="${2-}"; shift 2 ;;
    --dataset)       DATASET="${2-}"; shift 2 ;;
    --demand)        DEMANDS="${2-}"; DEMAND_PROVIDED=1; shift 2 ;;
    --demands)       DEMANDS="${2-}"; DEMAND_PROVIDED=1; shift 2 ;;
    --out)           OUT_CSV="${2-}"; shift 2 ;;
    -h|--help)
      echo "Usage:"
      echo "  $0 --all --dataset <dir-or-file> [--out results.csv]"
      echo "  $0 --solvers \"frt,ckr\" --dataset <dir-or-file> [--out results.csv]"
      echo "  $0 --solvers \"frt,ckr\" --dataset <dir-or-file> --demands \"gravity,gaussian\" [--out results.csv]"
      echo ""
      echo "All solvers and demands are passed in one Docker call per graph."
      echo "CSV has one row per (graph × solver × demand_model)."
      exit 0
      ;;
    *) echo "Unknown argument: $1"; exit 1 ;;
  esac
done

# Normalize repo-relative path mistake
if [[ "${DATASET:-}" == /* ]] && [[ "${DATASET:-}" == /experiments/* ]]; then
  DATASET="${DATASET#/}"
fi

if [[ "$RUN_ALL" -eq 1 ]]; then
  SOLVERS="$ALL_SOLVERS"
  DEMANDS="$ALL_DEMANDS"
  DEMAND_PROVIDED=1
fi

if [[ -z "${SOLVERS:-}" || -z "${DATASET:-}" ]]; then
  echo "ERROR: --solvers and --dataset are required (or use --all + --dataset)."
  exit 1
fi

TIMEOUT_BIN="${TIMEOUT_BIN:-timeout}"
if ! command -v "$TIMEOUT_BIN" >/dev/null 2>&1; then
  if command -v gtimeout >/dev/null 2>&1; then
    TIMEOUT_BIN="gtimeout"
  else
    echo "ERROR: 'timeout' not found. macOS: brew install coreutils"; exit 1
  fi
fi

mkdir -p "$(dirname "$OUT_CSV")"
mkdir -p "$OUT_DIR"

CSV="$OUT_CSV"

# Write header only if the file does not exist yet
if [[ ! -f "$CSV" ]]; then
  echo "dataset,graph,solver,num_nodes,num_edges,total_time_ms,solve_time_ms,transformation_time_ms,mwu_iterations,avg_oracle_time_ms,demand_model,offline_opt,achieved_congestion,ratio_pct,status" > "$CSV"
fi

# Collect graphs
DATASET_PATH="$DATASET"
GRAPHS=()
if [[ -f "$DATASET_PATH" ]]; then
  GRAPHS=("$DATASET_PATH")
  DATASET_ROOT="$(cd "$(dirname "$DATASET_PATH")" && pwd)"
else
  if [[ ! -d "$DATASET_PATH" ]]; then
    echo "ERROR: --dataset path not found: $DATASET_PATH"; exit 1
  fi
  DATASET_ROOT="$(cd "$DATASET_PATH" && pwd)"
  mapfile -t GRAPHS < <(find "$DATASET_ROOT" -type f -name "*.lgf" | sort)
fi

[[ ${#GRAPHS[@]} -gt 0 ]] || { echo "ERROR: No .lgf files found under $DATASET"; exit 1; }

# Normalise solvers/demands into a single comma-separated string each (no spaces)
SOLVERS_ARG="$(echo "$SOLVERS" | tr -d '[:space:]')"
DEMANDS_ARG="$(echo "${DEMANDS:-}" | tr -d '[:space:]')"

RUN_ID="$(date +%Y%m%d_%H%M%S)"

############################################
# Main loop — one Docker invocation per graph
############################################
for g in "${GRAPHS[@]}"; do
  g_abs="$(cd "$(dirname "$g")" && pwd)/$(basename "$g")"

  if [[ "$g_abs" == "$DATASET_ROOT/"* ]]; then
    rel_path="${g_abs#$DATASET_ROOT/}"
  else
    rel_path="$(basename "$g_abs")"
  fi

  if [[ "$rel_path" == *"/"* ]]; then
    dataset_label="${rel_path%%/*}"
  else
    dataset_label="$(basename "$DATASET_ROOT")"
  fi

  base="$(basename "$g_abs")"
  safe_rel="${rel_path//\//__}"
  log="$OUT_DIR/${safe_rel%.lgf}_${RUN_ID}.log"

  # Build argv for the container entrypoint: <solvers> <graph> [<demands>]
  cmd=( "$SOLVERS_ARG" "graphs/$rel_path" )
  if [[ "$DEMAND_PROVIDED" -eq 1 && -n "$DEMANDS_ARG" ]]; then
    cmd+=( "$DEMANDS_ARG" )
  fi

  echo "[RUN] $rel_path | solvers=$SOLVERS_ARG | demands=${DEMANDS_ARG:-none} | timeout=$TIMEOUT"

  status="OK"
  if "$TIMEOUT_BIN" --signal=SIGTERM --kill-after=30s "$TIMEOUT" \
      docker run --rm --runtime=runc \
        -e OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}" \
        -v "$DATASET_ROOT":/app/graphs \
        "$IMAGE" \
        "${cmd[@]}" \
      > "$log" 2>&1
  then
    status="OK"
  else
    ec=$?
    if   [[ "$ec" -eq 124 || "$ec" -eq 137 || "$ec" -eq 143 ]]; then status="TIMEOUT"
    elif [[ "$ec" -eq 139 ]]; then status="SEGFAULT"
    else status="ERROR_$ec"
    fi
  fi

  # ---------------------------------------------------------------------------
  # Parse the log — one row per (solver × demand_model)
  # ---------------------------------------------------------------------------
  awk \
    -v dataset="$dataset_label" \
    -v graph="$rel_path" \
    -v status="$status" \
    -v demands_arg="$DEMANDS_ARG" \
    -v demand_provided="$DEMAND_PROVIDED" \
  '
  BEGIN {
    nodes="NaN"; edges="NaN";
    solver=""; total_time="NaN"; solve_time="NaN"; transf_time="NaN";
    mwu="NaN"; avg_oracle="NaN";
    n_rows=0;
  }

  /^Graph loaded: [0-9]+ nodes, [0-9]+ edges\./ {
    tmp=$0
    sub(/^Graph loaded: /, "", tmp)
    nodes=tmp; sub(/ nodes.*$/, "", nodes)
    sub(/^[0-9]+ nodes, /, "", tmp)
    edges=tmp; sub(/ edges.*$/, "", edges)
    next
  }

  /^=== Running solver: / {
    solver=$0; sub(/^=== Running solver: /,"",solver); sub(/ ===/,"",solver)
    total_time="NaN"; solve_time="NaN"; transf_time="NaN";
    mwu="NaN"; avg_oracle="NaN";
    next
  }

  /^Total running time: [0-9.]+ ms$/ {
    tmp=$0; sub(/^Total running time: /,"",tmp); sub(/ ms$/,"",tmp); total_time=tmp; next
  }
  /^Solve time: [0-9.]+ ms$/ {
    tmp=$0; sub(/^Solve time: /,"",tmp); sub(/ ms$/,"",tmp); solve_time=tmp; next
  }
  /^Transformation time: [0-9.]+ ms$/ {
    tmp=$0; sub(/^Transformation time: /,"",tmp); sub(/ ms$/,"",tmp); transf_time=tmp; next
  }
  /^MWU iterations: [0-9]+$/ {
    tmp=$0; sub(/^MWU iterations: /,"",tmp); mwu=tmp; next
  }
  /^Average oracle time: [0-9.]+ ms$/ {
    tmp=$0; sub(/^Average oracle time: /,"",tmp); sub(/ ms$/,"",tmp); avg_oracle=tmp; next
  }

  /^Ratio off the optimal offline solution \[/ {
    dm=$0; sub(/^Ratio off the optimal offline solution \[/,"",dm); sub(/\] demand model:.*$/,"",dm)
    ratio_pct=$0; sub(/^.*: /,"",ratio_pct); sub(/%.*$/,"",ratio_pct)
    vals=$0; sub(/^.*\(/,"",vals); sub(/\).*$/,"",vals)
    n=split(vals, ab, " / ")
    offline_val  = (n>=1) ? ab[1] : "NaN"
    achieved_val = (n>=2) ? ab[2] : "NaN"
    printf "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n",
      dataset, graph, solver, nodes, edges,
      total_time, solve_time, transf_time, mwu, avg_oracle,
      dm, offline_val, achieved_val, ratio_pct, status
    n_rows++
    next
  }

  END {
    if (status != "OK" || n_rows == 0) {
      if (demand_provided == "1" && demands_arg != "") {
        n_d = split(demands_arg, dm_arr, ",")
      } else {
        n_d = 1; dm_arr[1] = "none"
      }
      if (n_rows == 0) {
        for (di=1; di<=n_d; di++) {
          printf "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n",
            dataset, graph,
            (solver=="" ? "unknown" : solver),
            nodes, edges,
            "NaN","NaN","NaN","NaN","NaN",
            dm_arr[di], "NaN","NaN","NaN", status
        }
      }
    }
  }
  ' "$log" >> "$CSV"

  echo "[DONE] $base | solvers=$SOLVERS_ARG | ${DEMANDS_ARG:-none} | $status"
done

echo "CSV written to: $CSV"
echo "Logs written to: $OUT_DIR"

