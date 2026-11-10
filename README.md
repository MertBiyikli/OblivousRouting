# E-Routing

### E-Routing is network resilience and routing analysis engine for evaluating robust routing strategies in complex networks. It provides tools for simulating network failures, analyzing routing protocols, and optimizing network performance under various conditions.

### Features
- **Solver ecosystem**: Multiple oblivious and semi-oblivious routing solvers
- **Traffic-demand evaluation**: Evaluate routing strategies under different traffic demands
- **Network failure simulation**: Simulate various network failure scenarios to assess routing robustness
- **Output formats**: JSON, text and console output Visualization export


### Build
Requirements:
- CMake >= 3.25
- C++23
- Boost libraries (optional, for certain features)
- Eigen3
- OR-Tools
- AMGCL


### Recommended build steps

Recommended build:
```
cmake --preset release
cmake --build --preset release
```
For tests:
```
cmake --preset debug-tests
cmake --build --preset debug-tests
ctest --preset debug --output-on-failure
```


### CLI

#### General usage:
```
./oblivious_routing solve \
--solver <solver> \
--graph <graph-file> \
[options]
```

Example:

```
./oblivious_routing solve \
--solver electrical \
--graph experiments/datasets/small/Backbone/1221.lgf
```
With demand evaluation:
```
./oblivious_routing solve \
--solver electrical \
--graph experiments/datasets/small/Backbone/1221.lgf \
--demand gravity
```

Compare multiple solvers:

```
./oblivious_routing solve \
--solver electrical,raecke_ckr,raecke_mst \
--graph experiments/datasets/small/Backbone/1221.lgf \
--demand uniform
```

Write results to JSON:

```
./oblivious_routing solve \
--solver electrical \
--graph experiments/datasets/small/Backbone/1221.lgf \
--output results/run.json
```


Now this info for developing only, for running the E-Routing studio, paste the following in the command line brather:

```
source .venv/bin/activate

export ER_NETWORK_REGISTRY="$PWD/e_routing/api/networks.json"
export ER_ENGINE="$PWD/build/mac-debug-tests-no-openmp/oblivious_routing"
export ER_IMPORTED_NETWORKS_DIR="$PWD/data/imported_networks"

python3 -m uvicorn e_routing.api.main:app \
  --host 127.0.0.1 \
  --port 8000 \
  --reload
```

#### Maintainer:
Mert Biyikli