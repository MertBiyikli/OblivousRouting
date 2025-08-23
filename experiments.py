#!/usr/bin/env python3
import argparse, subprocess, time, csv, sys, shlex
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from glob import glob
from datetime import datetime
from pathlib import Path
def guess_exe():
    names = [
        "cmake-build-debug/oblivious_routing",
        "cmake-build-relwithdebinfo/oblivious_routing",
        "cmake-build-release/oblivious_routing",
        "build/Debug/oblivious_routing",             # MSVC multi-config
        "build/Release/oblivious_routing",
    ]
    suf = ".exe" if (Path().anchor and Path().drive) else ""  # crude Windows check
    for n in names:
        p = Path(n + suf)
        if p.exists() and p.is_file():
            return str(p.resolve())
    return "./oblivious_routing" + suf

CANON = {
    # electrical
    "electrical":"electrical","electrical_naive":"electrical","ef":"electrical","e":"electrical","0":"electrical",
    # tree
    "tree":"tree","raecke":"tree","frt":"tree","r":"tree","t":"tree","1":"tree",
    # cohen
    "cohen":"cohen","lp":"cohen","applegate":"cohen","ac":"cohen","l":"cohen","2":"cohen",
}

def canon_solver(s: str) -> str:
    k = s.lower()
    if k not in CANON: raise ValueError(f"Unknown solver token: {s}")
    return CANON[k]

def run_one(exe: Path, solver: str, graph: Path, rep: int, timeout: float, logs_dir: Path):
    start = time.perf_counter()
    ts = datetime.now().isoformat(timespec="seconds")
    base = f"{solver}_{graph.stem}_rep{rep}"
    out_path = logs_dir / f"{base}.out"
    err_path = logs_dir / f"{base}.err"
    cmd = [str(exe), solver, str(graph)]
    try:
        with out_path.open("wb") as out, err_path.open("wb") as err:
            p = subprocess.run(cmd, stdout=out, stderr=err, timeout=timeout)
        timed_out = False
        rc = p.returncode
    except subprocess.TimeoutExpired:
        timed_out = True
        rc = 124
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
    }

def main():
    ap = argparse.ArgumentParser(description="Batch runner for oblivious_routing executable.")
    # in argparse:
    ap.add_argument("--exe", default=guess_exe(), help="Path to executable")
    ap.add_argument("--solvers", nargs="+", default=["electrical","tree","cohen"],
                    help="Solver tokens (e.g., electrical tree cohen or 0 1 2)")
    ap.add_argument("--graphs", nargs="+", required=True,
                    help="Graph files or globs (e.g., data/*.lfg data/small.lfg)")
    ap.add_argument("--repeats", type=int, default=1, help="Repeat each combo N times")
    ap.add_argument("--timeout", type=float, default=1800, help="Per-run timeout (seconds)")
    ap.add_argument("--jobs", type=int, default=1, help="Parallel workers")
    ap.add_argument("--out", default="results.csv", help="CSV output path")
    ap.add_argument("--logs", default="logs", help="Directory for stdout/stderr")
    args = ap.parse_args()

    exe = Path(args.exe).resolve()
    if not exe.exists():
        print(f"ERROR: executable not found: {exe}", file=sys.stderr)
        sys.exit(1)

    # Expand globs and sanity-check
    graphs = []
    for g in args.graphs:
        matches = sorted(Path(p) for p in glob(g))
        if not matches:
            print(f"WARNING: no files match {g}", file=sys.stderr)
        graphs.extend(matches)
    if not graphs:
        print("ERROR: no graph files to run.", file=sys.stderr)
        sys.exit(1)

    solvers = [canon_solver(s) for s in args.solvers]
    logs_dir = Path(args.logs)
    logs_dir.mkdir(parents=True, exist_ok=True)

    combos = []
    for s in solvers:
        for g in graphs:
            for r in range(args.repeats):
                combos.append((s, g, r))

    print(f"Running {len(combos)} jobs: {len(solvers)} solvers × {len(graphs)} graphs × {args.repeats} repeats")
    results = []
    with ThreadPoolExecutor(max_workers=max(1, args.jobs)) as pool:
        futs = [pool.submit(run_one, exe, s, g, r, args.timeout, logs_dir) for (s,g,r) in combos]
        for fut in as_completed(futs):
            res = fut.result()
            results.append(res)
            status = "TIMEOUT" if res["timed_out"] else f"rc={res['returncode']}"
            print(f"[{res['timestamp']}] {res['solver']} {Path(res['graph']).name}: {status}, {res['wall_time_s']}s")

    # Write CSV
    fieldnames = ["timestamp","solver","graph","returncode","timed_out","wall_time_s","stdout","stderr","cmd"]
    out_csv = Path(args.out)
    with out_csv.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in results:
            w.writerow(r)
    print(f"\nWrote {out_csv} with {len(results)} rows. Logs in {logs_dir}/")

if __name__ == "__main__":
    main()
