#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

# ======================
# CONFIG
# ======================

RESULTS_CSV = "results/results.csv"
OUT_DIR = Path("plots")
OUT_DIR.mkdir(exist_ok=True)

plt.rcParams.update({
    "font.size": 12,
    "axes.labelsize": 14,
    "legend.fontsize": 11,
    "xtick.labelsize": 11,
    "ytick.labelsize": 11,
    "figure.figsize": (6.5, 4.5),
})


# ======================
# LOAD DATA
# ======================

df = pd.read_csv(RESULTS_CSV)

required = {"solver", "num_edges", "total_time_ms", "mwu_iterations",
            "achieved_congestion", "offline_opt_value"}
assert required.issubset(df.columns), "Missing required columns in results.csv"


# Consistent solver order and colors
SOLVER_ORDER = sorted(df["solver"].unique())
COLORS = dict(zip(SOLVER_ORDER, plt.cm.tab10.colors))


# ======================
# A) Total running time vs number of edges
# ======================

plt.figure()

for solver in SOLVER_ORDER:
    sub = df[df["solver"] == solver].sort_values("num_edges")
    plt.plot(
        sub["num_edges"],
        sub["total_time_ms"],
        marker="o",
        label=solver,
        color=COLORS[solver],
    )

plt.xscale("log")
plt.yscale("log")
plt.xlabel("Number of edges")
plt.ylabel("Total running time [ms]")
plt.legend()
plt.tight_layout()
plt.savefig(OUT_DIR / "runtime_vs_edges.pdf")
plt.close()


# ======================
# B) MWU iterations vs number of edges
# ======================

plt.figure()

for solver in SOLVER_ORDER:
    sub = df[df["solver"] == solver].sort_values("num_edges")
    plt.plot(
        sub["num_edges"],
        sub["mwu_iterations"],
        marker="o",
        label=solver,
        color=COLORS[solver],
    )

plt.xscale("log")
plt.xlabel("Number of edges")
plt.ylabel("MWU iterations")
plt.legend()
plt.tight_layout()
plt.savefig(OUT_DIR / "mwu_iterations_vs_edges.pdf")
plt.close()


# ======================
# C) Error vs offline optimal
# ======================

df["relative_error"] = df["achieved_congestion"] / df["offline_opt_value"]

plt.figure()

for solver in SOLVER_ORDER:
    sub = df[df["solver"] == solver].sort_values("num_edges")
    plt.plot(
        sub["num_edges"],
        sub["relative_error"],
        marker="o",
        label=solver,
        color=COLORS[solver],
    )

plt.xscale("log")
plt.xlabel("Number of edges")
plt.ylabel("Relative error (solver / optimal)")
plt.legend()
plt.tight_layout()
plt.savefig(OUT_DIR / "error_vs_edges.pdf")
plt.close()

print("✔ Research-ready plots written to", OUT_DIR)
