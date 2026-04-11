#!/bin/bash

# Generate graphs with edge sizes up to 100,000 for each topology family

SCRIPT_DIR="/Users/halilibrahim/Desktop/Thesis/ObliviousRouting/scripts"
OUTPUT_BASE="/Users/halilibrahim/Desktop/Thesis/ObliviousRouting/experiments/datasets/synthetic"

mkdir -p "$OUTPUT_BASE/2dgrids"
mkdir -p "$OUTPUT_BASE/fattree"
mkdir -p "$OUTPUT_BASE/fatclique"
mkdir -p "$OUTPUT_BASE/expander"

echo "Generating 2D Grid topologies..."
# 2D Grid: edges ≈ 2*rows*cols - rows - cols
# Target edges: 10k, 20k, 30k, 40k, 50k, 60k, 70k, 80k, 90k, 100k
for i in 1 2 3 4 5 6 7 8 9 10; do
    rows=$((50 * i))
    cols=$((50 * i))
    cap=1
    python3 "$SCRIPT_DIR/generate_topologies.py" -t grid2d -o "$OUTPUT_BASE/2dgrids/grid2d_${i}.lgf" \
        --rows $rows --cols $cols --cap $cap --servers 0
done

echo "Generating Fat-Tree topologies..."
# Fat-Tree: vary switchRadix and numAgg
# Approximate edge formula: numAgg * (2 * (switchRadix/2)^2 + (switchRadix/2)^2 * switchRadix/2)
for i in 1 2 3 4 5 6 7 8 9 10; do
    switchRadix=$((4 + i * 2))
    numAgg=$((2 + i))
    python3 "$SCRIPT_DIR/generate_topologies.py" -t fattree_partial -o "$OUTPUT_BASE/fattree/fattree_${i}.lgf" \
        -p "{\"switchRadix\": $switchRadix, \"numAgg\": $numAgg}"
done

echo "Generating Fat-Clique topologies..."
# Fat-Clique: complete connections grow quadratically
# edges ≈ numLocalToR^2 * numSubblock^2 + numLocalToR^2 * numBlock^2 + numSubblock^2 * numBlock^2
for i in 1 2 3 4 5 6 7 8 9 10; do
    numServerPerToR=4
    numLocalToR=$((2 + i))
    numSubblock=$((2 + i))
    numBlock=$((2 + i))
    python3 "$SCRIPT_DIR/generate_topologies.py" -t fatclique -o "$OUTPUT_BASE/fatclique/fatclique_${i}.lgf" \
        -p "{\"numServerPerToR\": $numServerPerToR, \"numLocalToR\": $numLocalToR, \"numSubblock\": $numSubblock, \"numBlock\": $numBlock}"
done

echo "Generating Random Regular Expander topologies..."
# d-regular graph: edges = numNodes * degree / 2
# For target edges e, numNodes * degree = 2*e
# Let's vary: target 10k, 20k, ..., 100k edges
for i in 1 2 3 4 5 6 7 8 9 10; do
    target_edges=$((i * 10000))
    # Try different combinations of numNodes and degree
    # degree should be even for random_regular_graph
    degree=20
    numNodes=$((target_edges * 2 / degree))
    # Ensure numNodes * degree is even
    if [ $((numNodes * degree % 2)) -ne 0 ]; then
        numNodes=$((numNodes + 1))
    fi
    python3 "$SCRIPT_DIR/generate_topologies.py" -t expander -o "$OUTPUT_BASE/expander/expander_${i}.lgf" \
        -p "{\"numNodes\": $numNodes, \"degree\": $degree, \"linkCapacity\": 1, \"seed\": $i}"
done

echo "All graphs generated successfully!"

