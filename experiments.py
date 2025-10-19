#!/usr/bin/env python3
import argparse, subprocess, time, csv, sys, shlex
import concurrent
import os
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from glob import glob
from datetime import datetime
def count_edges(graph_path: Path) -> int:
    with graph_path.open() as f:
        lines = f.readlines()
    try:
        start = lines.index("@arcs\n")
    except ValueError:
        return 0
    count = 0
    for line in lines[start+1:]:
        if line.startswith("@"):  # next section begins
            break
        if line.strip():
            count += 1
    return count

def guess_exe():
    names = [
        #"cmake-build-debug/oblivious_routing",
        #"cmake-build-relwithdebinfo/oblivious_routing",
        "cmake-build-release/oblivious_routing", # remove parallel for sequential build
        "build/Debug/oblivious_routing",
        "build/Release/oblivious_routing",
    ]
    suf = ".exe" if (Path().anchor and Path().drive) else ""
    for n in names:
        p = Path(n + suf)
        if p.exists() and p.is_file():
            return str(p.resolve())
    return "./oblivious_routing" + suf

CANON = {
    "electrical": "electrical", "electrical_naive": "electrical", "ef": "electrical", "e": "electrical", "0": "electrical",
    "tree": "tree", "raecke": "tree", "frt": "tree", "r": "tree", "t": "tree", "1": "tree",
    "cohen": "cohen", "lp": "cohen", "applegate": "cohen", "ac": "cohen", "l": "cohen", "2": "cohen",
    "mst": "mst", "random_mst": "mst", "raecke_mst": "mst", "m": "mst", "3": "mst",
    "electrical_optimized": "electrical_optimized", "ef_opt": "electrical_optimized", "eo": "electrical_optimized", "4": "electrical_optimized",
}

CANON_DEMAND = {
    "gravity": "gravity", "gravity_model": "gravity",
    "uniform": "uniform", "uniform_model": "uniform",
    "gaussian": "gaussian", "gaussian_model": "gaussian",
    "bimodal": "bimodal", "bimodal_model": "bimodal",
}

def canon_solver(s: str) -> str:
    k = s.lower()
    if k not in CANON:
        raise ValueError(f"Unknown solver token: {s}")
    return CANON[k]

def canon_demand(s: str) -> str:
    k = s.lower()
    if k not in CANON_DEMAND:
        raise ValueError(f"Unknown demand token: {s}")
    return CANON_DEMAND[k]

import threading

# global: track processes per (solver, graph, rep)
group_processes = {}
group_lock = threading.Lock()

def run_one(exe: Path, solver: str, graph: Path, rep: int, timeout: float,
            logs_dir: Path, demand_model: str, number_edges: int = None) -> dict:
    start = time.perf_counter()
    ts = datetime.now().isoformat(timespec="seconds")
    # Construct output base name safely
    suffix = f"_{demand_model}" if demand_model else ""
    base = f"{solver}_{graph.stem}_rep{rep}{suffix}_{number_edges}_edges"
    out_path = logs_dir / f"{base}.out"
    err_path = logs_dir / f"{base}.err"

    # Build the command dynamically
    cmd = [str(exe), solver, str(graph)]
    if demand_model:
        cmd.append(demand_model)


    p = subprocess.Popen(cmd, stdout=out_path.open("wb"), stderr=err_path.open("wb"))

    # register this process under its group
    group_key = (solver, str(graph), rep)
    with group_lock:
        group_processes.setdefault(group_key, []).append(p)

    try:
        p.wait(timeout=timeout)
        timed_out = False
        rc = p.returncode
    except subprocess.TimeoutExpired:
        timed_out = True
        rc = 124
        # Kill all jobs in this group
        with group_lock:
            for proc in group_processes.get(group_key, []):
                if proc.poll() is None:  # still running
                    proc.kill()

    wall = time.perf_counter() - start
    return {
        "timestamp": ts,
        "solver": solver,
        "graph": str(graph),
        "returncode": rc,
        "timed_out": timed_out,
        "wall_time_s": f"{wall:.6f}",
        "stdout": str(out_path),
        "stderr": str(err_path),
        "cmd": " ".join(shlex.quote(x) for x in cmd),
        "demand": demand_model
    }


def main():
    ap = argparse.ArgumentParser(description="Batch runner for oblivious_routing executable.")
    ap.add_argument("--exe", default=guess_exe(), help="Path to executable")
    ap.add_argument("--solvers", nargs="+", default=["electrical", "tree", "cohen", "mst"],
                    help="Solver tokens (e.g., electrical tree cohen or 0 1 2 3)")
    ap.add_argument("--graphs", nargs="+", required=True,
                    help="Graph files or globs (e.g., data/*.lfg data/small.lfg)")
    ap.add_argument("--repeats", type=int, default=1, help="Repeat each combo N times")
    ap.add_argument("--timeout", type=float, default=1200, help="Per-run timeout (seconds)")
    ap.add_argument("--jobs", type=int, default=1, help="Parallel workers")
    ap.add_argument("--demand-model", nargs="+",  default=None, help="Demand model to use (e.g., gravity, binomial)")
    ap.add_argument("--dry-run", action="store_true", help="Only print jobs to run, don’t execute them")
    ap.add_argument("--smoke-test", action="store_true", help="Run only the first job to validate pipeline")

    ap.add_argument("-v", "--verbose", action="store_true", help="Verbose listing of jobs")
    output = datetime.now().strftime("result_%Y%m%d.csv")
    ap.add_argument("--out", default=output, help="CSV output path")
    ap.add_argument("--logs", default="logs", help="Directory for stdout/stderr")
    args = ap.parse_args()

    exe = Path(args.exe).resolve()
    if not exe.exists():
        print(f"ERROR: executable not found: {exe}", file=sys.stderr)
        sys.exit(1)

    graphs = []
    for g in args.graphs:
        matches = sorted(Path(p) for p in glob(g, recursive=True))
        if not matches:
            print(f"WARNING: no files match {g}", file=sys.stderr)
        graphs.extend(matches)
    if not graphs:
        print("ERROR: no graph files to run.", file=sys.stderr)
        sys.exit(1)

    solvers = [canon_solver(s) for s in args.solvers]

    if args.demand_model is not None:
        demands = [canon_demand(s) for s in args.demand_model]
    else:
        # make it empty to use default in executable
        demands = [None]

    logs_dir = Path(args.logs)
    logs_dir.mkdir(parents=True, exist_ok=True)

    combos = []
    for s in solvers:
        for g in graphs:
            for r in range(args.repeats):
                for d in demands:
                    combos.append((s, g, r, d))

    print(f"Running {len(combos)} jobs: {len(solvers)} solvers × {len(graphs)} graphs × {len(demands)} demand models × {args.repeats} repeats")
    if args.verbose:
        for i, (s, g, r, d) in enumerate(combos, start=1):
            print(f"  job {i:4}: solver={s} graph={g.name} rep={r} demand={d}")



    jobs = [
        (solver, graph, rep, args.timeout, logs_dir, demand)
        for solver in solvers
        for graph in graphs
        for rep in range(args.repeats)
        for demand in demands
    ]

    if args.dry_run:
        print(f"Running {len(jobs)} jobs: {len(solvers)} solvers × {len(graphs)} graphs × {args.repeats} repeats")
        for i, job in enumerate(jobs[:10]):  # only show a few
            solver, graph, rep, demand = job[:3]
            print(f"  job {i+1:4}: solver={solver} graph={os.path.basename(graph)} rep={rep} demand={demand}")

        if args.smoke_test and jobs:
            print("⚠️  Smoke test enabled: Running only the first job for verification...\n")
            result = run_one(exe, *jobs[0])
            print("Smoke test result:")
            print(result)
        sys.exit(0)

    all_ok = True
    results = []
    with ThreadPoolExecutor(max_workers=max(1, args.jobs)) as pool:

        futs = [pool.submit(run_one, exe, s, g, r, args.timeout, logs_dir, demand, count_edges(g))
                for (s, g, r, demand) in combos]

        for fut in as_completed(futs):
            try:
                res = fut.result()
                if res["returncode"] != 0 and all_ok:
                    all_ok = False
                    print(f"[!] A job failed with return code {res}")
                else:
                    results.append(res)
                status = "TIMEOUT" if res["timed_out"] else f"rc={res['returncode']}"
                print(f"[{res['timestamp']}] {res['solver']} {Path(res['graph']).name}: {status}, {res['wall_time_s']} demand: {res['demand']}")
            except Exception as e:
                print(f"[!] A job raised an exception: {e}", file=sys.stderr)
                all_ok = False

    # After all are done
    if all_ok:
        print("OK: All experiments completed successfully.")
    else:
        print("Some experiments failed.")


    # Write CSV
    fieldnames = ["timestamp", "solver", "graph", "returncode", "timed_out", "wall_time_s", "stdout", "stderr", "cmd", "demand"]
    out_csv = Path(args.out)
    with out_csv.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in results:
            w.writerow(r)
    print(f"\nWrote {out_csv} with {len(results)} rows. Logs in {logs_dir}/")

if __name__ == "__main__":
    main()
