#!/usr/bin/env bash
set -euo pipefail

OUT_CSV="results/results.csv"
OUT_DIR="${OUT_DIR:-results/logs}"
TIMEOUT="${TIMEOUT:-120m}"

# Native binary (override via: BIN=... ./run_experiments_using_binary.sh ...)
BIN="${BIN:-./build/oblivious_routing}"

SOLVERS=""
DATASET=""
DEMANDS=""
DEMAND_PROVIDED=0
CYCLE_REMOVAL="${CYCLE_REMOVAL:-none}"

ALL_SOLVERS="${ALL_SOLVERS:-electrical_naive, electrical_sketching,frt,ckr,mst,frt_mendel,ckr_mendel,frt_pointer,ckr_pointer,mst_pointer,frt_mendel_pointer,ckr_mendel_pointer,lp}"
ALL_DEMANDS="${ALL_DEMANDS:-gravity,gaussian,uniform,bimodal}"

RUN_ALL=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --all|-all)      RUN_ALL=1; shift ;;
    --solvers)       SOLVERS="${2-}"; shift 2 ;;
    --dataset)       DATASET="${2-}"; shift 2 ;;
    --demand)        DEMANDS="${2-}"; DEMAND_PROVIDED=1; shift 2 ;;
    --demands)       DEMANDS="${2-}"; DEMAND_PROVIDED=1; shift 2 ;;
    --cycle-removal) CYCLE_REMOVAL="${2-}"; shift 2 ;;
    --bin)           BIN="${2-}"; shift 2 ;;
    --out)           OUT_CSV="${2-}"; shift 2 ;;
    -h|--help)
      echo "Usage:"
      echo "  $0 --all --dataset <dir-or-file> [--out results.csv]"
      echo "  $0 --solvers \"frt,ckr\" --dataset <dir-or-file> [--out results.csv]"
      echo "  $0 --solvers \"frt,ckr\" --dataset <dir-or-file> --demands \"gravity,gaussian\" [--out results.csv]"
      echo "  $0 --solvers \"frt,ckr\" --dataset <dir-or-file> --cycle-removal <strategy> [--out results.csv]"
      echo ""
      echo "Cycle removal strategies: none, tarjan, linear (default: none)"
      echo ""
      echo "All solvers, demands, and cycle removal strategy are passed in one binary call per graph."
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

if [[ ! -x "$BIN" ]]; then
  echo "ERROR: binary not found or not executable: $BIN"
  echo "Hint: cmake -S . -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build -j"
  exit 1
fi

mkdir -p "$(dirname "$OUT_CSV")"
mkdir -p "$OUT_DIR"

CSV="$OUT_CSV"

# Always write a fresh header — each invocation owns its output file
echo "dataset,graph,solver,num_nodes,num_edges,total_time_micro_seconds,solve_time_micro_seconds,transformation_time_micro_seconds,mwu_iterations,avg_oracle_time_micro_seconds,avg_tree_height,load_computation_micro_seconds,mendel_total_micro_seconds,mendel_avg_micro_seconds,oblivious_ratio,demand_model,offline_opt,achieved_congestion,ratio_pct,mwu_weight_update_time_micro_seconds,cycle_removal_type,cycle_removal_time_micro_seconds,status" > "$CSV"

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
# Main loop — one binary invocation per graph
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

  # Build command: binary <solvers> <graph> [<demands>] [<cycle-removal>]
  cmd=( "$BIN" "$SOLVERS_ARG" "$g_abs" )
  if [[ "$DEMAND_PROVIDED" -eq 1 && -n "$DEMANDS_ARG" ]]; then
    cmd+=( "$DEMANDS_ARG" )
  fi
  if [[ -n "$CYCLE_REMOVAL" && "$CYCLE_REMOVAL" != "none" ]]; then
    cmd+=( "$CYCLE_REMOVAL" )
  fi

  echo "[RUN] $rel_path | solvers=$SOLVERS_ARG | demands=${DEMANDS_ARG:-none} | cycle_removal=$CYCLE_REMOVAL | timeout=$TIMEOUT"

  status="OK"
  if "$TIMEOUT_BIN" --signal=SIGTERM --kill-after=30s "$TIMEOUT" \
      env OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}" \
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
  # Parse the log.
  #
  # Output structure (one binary call, multiple solvers, multiple demands):
  #
  #   Graph loaded: 36 nodes, 96 edges.
  #   ...
  #   === Running solver: raecke_frt ===
  #   Total running time: 95 ms
  #   Solve time: 50 ms
  #   Transformation time: 24 ms
  #   MWU iterations: 21
  #   Average oracle time: 2.28571 ms
  #   Ratio off the optimal offline solution [bimodal] demand model: 519.863% (2014.49 / 10472.6)
  #   Ratio off the optimal offline solution [uniform] demand model: 717.007% (11557 / 82864.5)
  #   ...
  #   === Running solver: raecke_ckr ===
  #   ...
  #
  # Strategy: track current solver section; for every ratio line emit one CSV row.
  # If status != OK, emit one row per (solver x demand) with NaN values.
  # ---------------------------------------------------------------------------
  # Parse the log into a temp file, then append atomically to the CSV.
  # Using a temp file ensures no partial/stale data from a previous run
  # can leak into the output even if the script is interrupted.
  # ---------------------------------------------------------------------------
  _tmp_rows="$(mktemp)"
  awk \
    -v dataset="$dataset_label" \
    -v graph="$rel_path" \
    -v status="$status" \
    -v demands_arg="$DEMANDS_ARG" \
    -v demand_provided="$DEMAND_PROVIDED" \
  '
  function flush_no_demand_row() {
    printf "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n",
      dataset, graph, solver, nodes, edges,
      total_time, solve_time, transf_time, mwu, avg_oracle,
      avg_tree_height, load_computation_micro_seconds, mendel_total, mendel_avg, oblivious_ratio,
      "none", "NaN","NaN","NaN", mwu_weight_update_time, cycle_removal_type, cycle_removal_time, status
    n_rows++
  }
  function reset_solver_state() {
    total_time="NaN"; solve_time="NaN"; transf_time="NaN";
    mwu="NaN"; avg_oracle="NaN"; avg_tree_height="NaN";
    mendel_total="NaN"; mendel_avg="NaN"; oblivious_ratio="NaN";
    mwu_weight_update_time="NaN";
    load_computation_micro_seconds="NaN";
    cycle_removal_type="none"; cycle_removal_time="NaN";
  }
  BEGIN {
    nodes="NaN"; edges="NaN";
    solver=""; n_rows=0; solver_seen=0;
    mwu_weight_update_time="NaN";
    load_computation_micro_seconds="NaN";
    cycle_removal_type="none"; cycle_removal_time="NaN";
    reset_solver_state()
  }

  # Graph metadata (appears once at the top)
  /^Graph loaded: [0-9]+ nodes, [0-9]+ edges\./ {
    tmp=$0
    sub(/^Graph loaded: /, "", tmp)
    nodes=tmp; sub(/ nodes.*$/, "", nodes)
    sub(/^[0-9]+ nodes, /, "", tmp)
    edges=tmp; sub(/ edges.*$/, "", edges)
    next
  }

  # Cycle cancellation strategy
  /^Using cycle cancellation strategy: / {
    tmp=$0; sub(/^Using cycle cancellation strategy: /, "", tmp); cycle_removal_type=tmp; next
  }

  # Cycle cancellation time
  /^Cycle cancellation time: [0-9][0-9.]*([eE][+-]?[0-9]+)?  *micro  *seconds/ {
    tmp=$0; sub(/^Cycle cancellation time: /, "", tmp); sub(/ *micro  *seconds/, "", tmp); cycle_removal_time=tmp; next
  }

  # New solver section — flush previous solver row if no demands, then reset
  /^=== Running solver: / {
    if (demand_provided != "1" && solver != "" && solver_seen) {
      flush_no_demand_row()
    }
    solver=$0; sub(/^=== Running solver: /,"",solver); sub(/ ===/,"",solver)
    reset_solver_state()
    solver_seen=1
    next
  }

  /^Total time: [0-9][0-9.]*([eE][+-]?[0-9]+)? micro[ _]seconds/ {
    tmp=$0; sub(/^Total time: /,"",tmp); sub(/ micro[ _]seconds/,"",tmp); total_time=tmp; next
  }
  /^Solve time: [0-9][0-9.]*([eE][+-]?[0-9]+)? micro[ _]seconds/ {
    tmp=$0; sub(/^Solve time: /,"",tmp); sub(/ micro[ _]seconds/,"",tmp); solve_time=tmp; next
  }
  /^Transformation time: [0-9][0-9.]*([eE][+-]?[0-9]+)? micro[ _]seconds/ {
    tmp=$0; sub(/^Transformation time: /,"",tmp); sub(/ micro[ _]seconds/,"",tmp); transf_time=tmp; next
  }
  /^MWU iterations: [0-9]+$/ {
    tmp=$0; sub(/^MWU iterations: /,"",tmp); mwu=tmp; next
  }
  /^MWU load computation: [0-9][0-9.]*([eE][+-]?[0-9]+)? micro[ _]seconds/ {
      tmp=$0; sub(/^MWU load computation: /,"",tmp); sub(/ micro[ _]seconds/,"",tmp); load_computation_micro_seconds=tmp; next
    }
  /^Average oracle time: [0-9][0-9.]*([eE][+-]?[0-9]+)? micro[ _]seconds/ {
    tmp=$0; sub(/^Average oracle time: /,"",tmp); sub(/ micro[ _]seconds/,"",tmp); avg_oracle=tmp; next
  }
  /^Total MWU weight update time: [0-9][0-9.]*([eE][+-]?[0-9]+)? micro[ _]seconds/ {
    tmp=$0; sub(/^Total MWU weight update time: /,"",tmp); sub(/ micro[ _]seconds/,"",tmp); mwu_weight_update_time=tmp; next
  }
  /^Average tree height: [0-9][0-9.]*([eE][+-]?[0-9]+)?$/ {
    tmp=$0; sub(/^Average tree height: /,"",tmp); avg_tree_height=tmp; next
  }
  /^Total time spent on Mendel scaling: [0-9][0-9.]*([eE][+-]?[0-9]+)?  *micro  *seconds/ {
    tmp=$0; sub(/^Total time spent on Mendel scaling: /,"",tmp); sub(/ *micro  *seconds/,"",tmp); mendel_total=tmp; next
  }
  /^Average time spent on Mendel scaling per iteration: [0-9][0-9.]*([eE][+-]?[0-9]+)?  *micro  *seconds/ {
    tmp=$0; sub(/^Average time spent on Mendel scaling per iteration: /,"",tmp); sub(/ *micro  *seconds/,"",tmp); mendel_avg=tmp; next
  }
  /^Oblivious ratio: [0-9][0-9.]*([eE][+-]?[0-9]+)?$/ {
    tmp=$0; sub(/^Oblivious ratio: /,"",tmp); oblivious_ratio=tmp; next
  }

  # Ratio line — one row per (solver × demand_model)
  /^Ratio off the optimal offline solution \[/ {
    dm=$0; sub(/^Ratio off the optimal offline solution \[/,"",dm); sub(/\] demand model:.*$/,"",dm)
    ratio_pct=$0; sub(/^.*: /,"",ratio_pct); sub(/%.*$/,"",ratio_pct)
    vals=$0; sub(/^.*\(/,"",vals); sub(/\).*$/,"",vals)
    n=split(vals, ab, " / ")
    offline_val  = (n>=1) ? ab[1] : "NaN"
    achieved_val = (n>=2) ? ab[2] : "NaN"
    printf "\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",\"%s\"\n",
      dataset, graph, solver, nodes, edges,
      total_time, solve_time, transf_time, mwu, avg_oracle,
      avg_tree_height, load_computation_micro_seconds, mendel_total, mendel_avg, oblivious_ratio,
      dm, offline_val, achieved_val, ratio_pct, mwu_weight_update_time, cycle_removal_type, cycle_removal_time, status
    n_rows++
    next
  }

  END {
    if (n_rows == 0 && !solver_seen) {
      # Binary produced no output at all (crash before any solver ran)
      if (demand_provided == "1" && demands_arg != "") {
        n_d = split(demands_arg, dm_arr, ",")
      } else {
        n_d = 1; dm_arr[1] = "none"
      }
      for (di=1; di<=n_d; di++) {
        printf "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n",
          dataset, graph, "unknown", nodes, edges,
          "NaN","NaN","NaN","NaN","NaN","NaN","NaN","NaN","NaN",
          dm_arr[di], "NaN","NaN","NaN", mwu_weight_update_time, cycle_removal_type, cycle_removal_time, status
      }
    } else if (demand_provided != "1" && solver_seen) {
      # Flush the last solver (no demand model run)
      flush_no_demand_row()
    }
  }
  ' "$log" > "$_tmp_rows"
  cat "$_tmp_rows" >> "$CSV"
  rm -f "$_tmp_rows"

  echo "[DONE] $base | solvers=$SOLVERS_ARG | ${DEMANDS_ARG:-none} | $CYCLE_REMOVAL | $status"
done

echo "CSV written to: $CSV"
echo "Logs written to: $OUT_DIR"

