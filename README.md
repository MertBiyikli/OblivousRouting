## Oblivious Routing algorithms

This repository contains implementations of various Oblivious Routing algorithms. Oblivious Routing is a technique used in network routing where the path taken by packets is determined without knowledge of the current network state or traffic conditions.
More precisely, we implemented an electrical flow-based and tree-based algorithm for solving the Oblivious Routing problem.
The observed algorithms include:

**MWU-based Oblivious Routing Approximation Algorithms:**
- Räcke STOC 2008 
  **Optimal Hierarchical Decompositions for Congestion Minimization in Networks**[[paper]](https://dl.acm.org/doi/10.1145/1374376.1374415)
- Goranci et al.  ITCS 2024
  **Electrical Flows for Polylogarithmic Competitive Oblivious Routing**[[paper]](https://arxiv.org/abs/2303.02491)

LP-based evaluation of the optimal Oblivious Routing:
- Applegate and Cohen SIGCOMM 2003 **Making intra-domain routing robust to changing and uncertain traffic demands: understanding fundamental tradeoffs**[[paper]](https://dl.acm.org/doi/10.1145/863955.863991)


Tree-based oracle variants implemented in this repository include:
- **FRT** — Fakcharoenphol, Rao, Talwar (2004) **A tight bound on approximating arbitrary metrics by tree metrics**[[paper]](https://dl.acm.org/doi/abs/10.1145/780542.780608)
- **Fast-CKR** — Mendel and Schwob (2009) **Fast C-K-R Partitions of Sparse Graphs ∗**[[paper]](https://arxiv.org/abs/0809.1902)
- **MST-based** decompositions

In addition we provide the following scaling refinement mechanism to the tree-based oracle variants:
- **MendelScaling** - Mendel and Schwob (2009): scales down the required number of scales when computing the HST while retaining the stretch quality

Electrical flow based oracle variant is based on the [[AMGCL solver]](https://github.com/ddemidov/amgcl).
The default solver is set to use the algebraic multigrid method for preconditioning the Conjugate Gradient (CG) method.

We provide two variants for the electrical flow solver:
 - **Computing the exact loads**: This variant computes the exact loads on the edges of the network based on the electrical flow model.
 - **Approximating the loads**: This variant uses an approximation method to estimate the loads on the edges while being theoretically sound and empirically effective. This approach can significantly reduce the computational overhead while still providing good performance in practice.

### Repository Structure

- `app/main.cpp`: Contains the main application code for running the Oblivious Routing algorithms and evaluating their performance.
- `include/`: Contains header files for the Oblivious Routing algorithms and related utilities.
- `source/`: Contains the source code for the Oblivious Routing algorithms.
- `experiments/datasets`: Contains datasets and network topologies used for testing the algorithms.


### Getting Started
#### Prerequisites
- C++ compiler (supporting C++11 or later)
- CMake (version 3.20 or later)
- Boost Libraries
- Eigen3
- OR-Tools (for LP solving)

#### Building the Project (Recommended)

This project provides predefined CMake presets for reproducible builds.

From the project root, simply run:

```bash
cmake --preset release
cmake --build --preset release
```

Change the preset configuration from `release` to `debug` for a debug build.

   
#### Running the Algorithms
- To run an oblivious Routing algorithms, use the following command:
   ```bash
   ./oblivious_routing <algorithm> <graph> [OPTIONAL]:<demand_model> <graph_format>
   ```
   Replace `<graph>` with the path to your network topology file and `<algorithm>` with `"electrical,frt,ckr,mst,cohen"`.
   Replace `<demand_model>` with the desired traffic demand model (e.g., `uniform, gravity, bimodal, gaussian`). The demand model 
   evaluates how traffic is distributed across the network given the precomputed routing scheme.


### Running with Docker
You can also run the project using Docker. First, build the Docker image:
```bash
    docker build -t oblivious-routing .
```

For a simplified execution, you can run the Docker container with the following script
```bash
    ./scripts/run_experiment.sh --solvers <algorithm> --dataset <graph> --demand <demand_model> --out <output_directory>
```
### Example smoke test
To perform a quick smoke test, you can run the following command:
```bash
    ./scripts/run_experiment.sh --solvers "electrical,frt,ckr,mst,cohen" --dataset /experiments/datasets/small/Backbone/1221.lgf --demand uniform --out results/test.csv
```


##### Contributors (sorted alphabetically by last name):
- Mert Biyikli (maintainer)