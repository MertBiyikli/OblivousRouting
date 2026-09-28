#!/usr/bin/env bash
set -euo pipefail

OUT_CSV="results/results.csv"
OUT_DIR="${OUT_DIR:-results/logs}"
TIMEOUT="${TIMEOUT:-120m}"

# Native binary (override via: BIN=... ./run_experiments_using_binary.sh ...)
BIN="${BIN:-../cmake-build-release/oblivious_routing}"

SOLVERS=""
DATASET=""
DEMANDS=""
DEMAND_PROVIDED=0
OUTPUT_FORMAT="${OUTPUT_FORMAT:-cout}"

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
    --bin)           BIN="${2-}"; shift 2 ;;
    --out)           OUT_CSV="${2-}"; shift 2 ;;
    --format)        OUTPUT_FORMAT="${2-}"; shift 2 ;;
    -h|--help)
      echo "Usage:"
      echo "  $0 --all --dataset <dir-or-file> [--out results.csv] [--format cout]"
      echo "  $0 --solvers \"frt,ckr\" --dataset <dir-or-file> [--out results.csv] [--format cout]"
      echo "  $0 --solvers \"frt,ckr\" --dataset <dir-or-file> --demands \"gravity,gaussian\" [--out results.csv] [--format cout]"
      echo "All solvers and demands are passed in one binary call per graph."
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

# Always write a fresh header — each invocation owns its output file.
# The parser below uses the same schema for every solver type; unavailable
# solver-specific metrics are written as NaN.
CSV_FIELDS="date,dataset,graph,solver,routing_base,num_nodes,num_edges,execution_status,solver_status,total_time_microseconds,preprocessing_time_microseconds,solve_time_microseconds,oblivious_ratio,demand_model,demand_congestion,demand_runtime_microseconds,offline_opt,ratio_pct,candidate_paths,average_paths_per_pair,mwu_iterations,mwu_solve_time_microseconds,mwu_transformation_time_microseconds,mwu_load_computation_time_microseconds,mwu_weight_update_time_microseconds,mwu_average_oracle_time_microseconds,mwu_oracle_calls,hierarchy_runtime_microseconds,tree_construction_runtime_microseconds,basis_flow_runtime_microseconds,hierarchy_levels,hierarchy_clusters,clusters_per_level,maximum_cluster_vertices,average_cluster_vertices,tree_nodes,tree_edges,tree_depth,basis_flows,total_electrical_solves,average_electrical_solves_per_basis_flow,maximum_basis_embedding_congestion,maximum_conservation_error,average_tree_height,mendel_total_microseconds,mendel_average_microseconds"
echo "$CSV_FIELDS" > "$CSV"

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
  graph_label="$(basename "$(dirname "$g_abs")")/$base"
  safe_rel="${rel_path//\//__}"
  log="$OUT_DIR/${safe_rel%.lgf}_${RUN_ID}.log"

  # Build command: binary <solvers> <graph> [<demands>] <output-format>
  cmd=( "$BIN" "$SOLVERS_ARG" "$g_abs" )
  if [[ "$DEMAND_PROVIDED" -eq 1 && -n "$DEMANDS_ARG" ]]; then
    cmd+=( "$DEMANDS_ARG" )
  fi
  cmd+=( "$OUTPUT_FORMAT" )


  echo "[RUN] $rel_path | solvers=$SOLVERS_ARG | demands=${DEMANDS_ARG:-none} | format=$OUTPUT_FORMAT | timeout=$TIMEOUT"

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

  # Parse both the current labelled cout format and the historical
  # "=== Running solver" / ratio-line format. Rows are buffered until the end
  # of each solver section because MWU and hierarchy metrics follow demands.
  _tmp_rows="$(mktemp)"
  python3 - "$log" "$dataset_label" "$graph_label" "$status" \
    "$DEMANDS_ARG" "$DEMAND_PROVIDED" "$CSV_FIELDS" > "$_tmp_rows" <<'PY'
import csv
import re
import sys

log_path, dataset, fallback_graph, execution_status, demands_arg, demand_provided, fields_arg = sys.argv[1:]
fields = fields_arg.split(",")
NA = "NaN"

label_to_field = {
    "Total runtime (microseconds)": "total_time_microseconds",
    "Preprocessing runtime (microseconds)": "preprocessing_time_microseconds",
    "Solve runtime (microseconds)": "solve_time_microseconds",
    "Oblivious ratio": "oblivious_ratio",
    "Candidate paths": "candidate_paths",
    "Average paths per pair": "average_paths_per_pair",
    "Iterations": "mwu_iterations",
    "Solve time (microseconds)": "mwu_solve_time_microseconds",
    "Transformation time (microseconds)": "mwu_transformation_time_microseconds",
    "Load computation time (microseconds)": "mwu_load_computation_time_microseconds",
    "Weight update time (microseconds)": "mwu_weight_update_time_microseconds",
    "Average oracle time (microseconds)": "mwu_average_oracle_time_microseconds",
    "Oracle calls": "mwu_oracle_calls",
    "Hierarchy runtime (microseconds)": "hierarchy_runtime_microseconds",
    "Tree construction runtime (microseconds)": "tree_construction_runtime_microseconds",
    "Basis-flow runtime (microseconds)": "basis_flow_runtime_microseconds",
    "Hierarchy levels": "hierarchy_levels",
    "Hierarchy clusters": "hierarchy_clusters",
    "Clusters per level": "clusters_per_level",
    "Maximum cluster vertices": "maximum_cluster_vertices",
    "Average cluster vertices": "average_cluster_vertices",
    "Tree nodes": "tree_nodes",
    "Tree edges": "tree_edges",
    "Tree depth": "tree_depth",
    "Basis flows": "basis_flows",
    "Total electrical solves": "total_electrical_solves",
    "Average electrical solves per basis flow": "average_electrical_solves_per_basis_flow",
    "Maximum basis embedding congestion": "maximum_basis_embedding_congestion",
    "Maximum conservation error": "maximum_conservation_error",
}

def fresh():
    row = {field: NA for field in fields}
    row.update(dataset=dataset, graph=fallback_graph,
               execution_status=execution_status)
    return {"row": row, "demands": [], "current_demand": None}

sections = []
section = None

def ensure_section():
    global section
    if section is None:
        section = fresh()
    return section

def finish_section():
    global section
    if section is not None and section["row"]["solver"] != NA:
        sections.append(section)
    section = None

with open(log_path, encoding="utf-8", errors="replace") as handle:
    for raw in handle:
        line = raw.strip()
        if not line:
            continue

        if line.startswith("Date: "):
            finish_section()
            ensure_section()["row"]["date"] = line[6:].strip()
            continue

        old_solver = re.fullmatch(r"=== Running solver: (.*?) ===", line)
        if old_solver:
            finish_section()
            ensure_section()["row"]["solver"] = old_solver.group(1)
            continue

        if line.startswith("Solver: "):
            ensure_section()["row"]["solver"] = line[8:].strip()
            continue

        sec = ensure_section()
        row = sec["row"]

        graph_loaded = re.fullmatch(r"Graph loaded: (\d+) nodes, (\d+) edges\.", line)
        if graph_loaded:
            row["num_nodes"], row["num_edges"] = graph_loaded.groups()
            continue

        demand = re.fullmatch(r"Demand \[(.+)]", line)
        if demand:
            entry = {"demand_model": demand.group(1),
                     "demand_congestion": NA,
                     "demand_runtime_microseconds": NA,
                     "offline_opt": NA, "ratio_pct": NA}
            sec["demands"].append(entry)
            sec["current_demand"] = entry
            continue

        if sec["current_demand"] is not None:
            if line.startswith("Congestion: "):
                sec["current_demand"]["demand_congestion"] = line.split(":", 1)[1].strip()
                continue
            if line.startswith("Runtime (microseconds): "):
                sec["current_demand"]["demand_runtime_microseconds"] = line.split(":", 1)[1].strip()
                continue

        ratio = re.fullmatch(
            r"Ratio off the optimal offline solution \[(.+?)] demand model: "
            r"([^%]+)% \(([^/]+) / ([^)]+)\)", line)
        if ratio:
            dm, ratio_pct, offline_opt, congestion = (value.strip() for value in ratio.groups())
            sec["demands"].append({"demand_model": dm,
                "demand_congestion": congestion, "demand_runtime_microseconds": NA,
                "offline_opt": offline_opt, "ratio_pct": ratio_pct})
            continue

        simple = {
            "Nodes": "num_nodes", "Edges": "num_edges",
            "Routing base": "routing_base", "Status": "solver_status",
            "Total time": "total_time_microseconds",
            "Solve time": "solve_time_microseconds",
            "Transformation time": "mwu_transformation_time_microseconds",
            "MWU iterations": "mwu_iterations",
            "MWU load computation": "mwu_load_computation_time_microseconds",
            "Average oracle time": "mwu_average_oracle_time_microseconds",
            "Total MWU weight update time": "mwu_weight_update_time_microseconds",
            "Average tree height": "average_tree_height",
            "Total time spent on Mendel scaling": "mendel_total_microseconds",
            "Average time spent on Mendel scaling per iteration": "mendel_average_microseconds",
        }
        matched = False
        for label, field in simple.items():
            if line.startswith(label + ": "):
                value = line[len(label) + 2:].strip()
                value = re.sub(r"\s+micro[ _ ]*seconds$", "", value)
                row[field] = value
                matched = True
                break
        if matched:
            continue

        for label, field in label_to_field.items():
            if line.startswith(label + ": "):
                row[field] = line[len(label) + 2:].strip()
                break

finish_section()

if not sections:
    failed = fresh()
    failed["row"]["solver"] = "unknown"
    requested = demands_arg.split(",") if demand_provided == "1" and demands_arg else ["none"]
    failed["demands"] = [{"demand_model": dm, "demand_congestion": NA,
                           "demand_runtime_microseconds": NA,
                           "offline_opt": NA, "ratio_pct": NA} for dm in requested]
    sections.append(failed)

writer = csv.DictWriter(sys.stdout, fieldnames=fields, lineterminator="\n")
for sec in sections:
    demands = sec["demands"] or [{"demand_model": "none"}]
    for demand in demands:
        output = sec["row"].copy()
        output.update(demand)
        writer.writerow(output)
PY
  cat "$_tmp_rows" >> "$CSV"
  rm -f "$_tmp_rows"

  echo "[DONE] $base | solvers=$SOLVERS_ARG | ${DEMANDS_ARG:-none} | $status"
done

echo "CSV written to: $CSV"
echo "Logs written to: $OUT_DIR"