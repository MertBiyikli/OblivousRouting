## Oblivious Routing algorithms

This repository contains implementations of various Oblivious Routing algorithms. Oblivious Routing is a technique used in network routing where the path taken by packets is determined without knowledge of the current network state or traffic conditions.
The goal is to compare two MWU style algorithm each instantiating different oracles.
More precisely, we examine how electrical flow and tree based oracles solve the Oblivious Routing.
The observed algorithms include:

MWU-based Oblivious Routing Approximation Algorithms:
- Räcke STOC 2008 ("Optimal Hierarchical Decompositions for Congestion Minimization in Networks")
- Goranci et al.  ITCS 2024 ("Electrical Flows for Polylogarithmic Competitive Oblivious Routing")

LP-based evaluation of the optimal Oblivious Routing:
- Applegate and Cohen SIGCOMM 2003 ("Making intra-domain routing robust to changing and uncertain traffic demands: understanding fundamental tradeoffs")


Note that we have implement different variants of the tree based algorithms including FRT, CKR, and Minimum Spanning Tree based oracles.

- FRT = Fakcharoenphol, Rao, and Talwar 2004 ("A tight bound on approximating arbitrary metrics by tree metrics")
- CKR = Mendel and Schwob 2009 ("Fast C-K-R Partitions of Sparse Graphs ∗")

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
- To run the MWU-based Oblivious Routing algorithms, use the following command:
   ```bash
   ./oblivious_routing <algorithm> <graph> <demand_model>
   ```
   Replace `<graph>` with the path to your network topology file and `<algorithm>` with `"electrical,frt,ckr,mst,cohen"`.
   Replace `<demand_model>` with the desired traffic demand model (e.g., `uniform, gravity, bimodal, gaussian`).