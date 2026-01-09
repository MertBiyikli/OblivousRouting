## Oblivious Routing algorithms

This repository contains implementations of various Oblivious Routing algorithms. Oblivious Routing is a technique used in network routing where the path taken by packets is determined without knowledge of the current network state or traffic conditions.
The goal is to compare two MWU style algorithm each instantiating different oracles.
More precisely, we examine how electrical flow and tree based oracles solve the Oblivious Routing.
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
- **CKR** — Mendel and Schwob (2009) **Fast C-K-R Partitions of Sparse Graphs ∗**[[paper]](https://arxiv.org/abs/0809.1902)
- **MST-based** decompositions


### Repository Structure

- `src/`: Contains the source code for the Oblivious Routing algorithms.
- `data/`: Contains datasets and network topologies used for testing the algorithms.


### Getting Started
#### Prerequisites
- C++ compiler (supporting C++11 or later)
- CMake (version 3.20 or later)
- Boost Libraries
- Eigen3
- OR-Tools (for LP solving)

#### Building the Project
1. Clone the repository:
   ```bash
   git clone "https://github.com/MertBiyikli/OblivousRouting.git"
   cd OblivousRouting
   mkdir build
   cd build
   cmake .. -DCMAKE_BUILD_TYPE=Release
   make -j
    ```
   
#### Running the Algorithms
- To run an oblivious Routing algorithms, use the following command:
   ```bash
   ./oblivious_routing <algorithm> <graph> <demand_model>
   ```
   Replace `<graph>` with the path to your network topology file and `<algorithm>` with `"electrical,frt,ckr,mst,cohen"`.
   Replace `<demand_model>` with the desired traffic demand model (e.g., `uniform, gravity, bimodal, gaussian`). The demand model 
   evaluates how traffic is distributed across the network givent the precomputed routing scheme.


### Running with Docker (recommended)
You can also run the project using Docker. First, build the Docker image:
```bash
    docker build -t oblivious_routing .
```

For a simplified execution, you can run the Docker container with the following script
```bash
    ./scripts/run_experiment.sh --solvers <algorithm> --dataset <graph> --demand <demand_model> --out <output_directory>
```
### Example smoke test
To perform a quick smoke test, you can run the following command:
```bash
    ./scripts/run_experiment.sh --solvers "electrical,frt,ckr,mst,cohen" --dataset /experiments/datasets/Backbone/1221.lgf --demand uniform --out results/test.csv
```


##### Contributors (sorted alphabetically by last name):
- Mert Biyikli (maintainer)
- Gramoz Goranci