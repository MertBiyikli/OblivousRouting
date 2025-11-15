//
// Created by Mert Biyikli on 25.10.25.
//

#ifndef OBLIVIOUSROUTING_RAECKE_CKR_SOLVER_H
#define OBLIVIOUSROUTING_RAECKE_CKR_SOLVER_H



#include "../raecke_base.h"
#include "raecke_ckr_transform.h"
#include "ckr_tree_decomposer.h"



class CKRSolver : public RaeckeBase<TreeNode*> {
public:
    void init(Graph& g) override {
        m_graph = g;
    }


    TreeNode* getTree(Graph& g) override{
        g.createDistanceMatrix();
        double delta = g.GetDiameter();
        TreeDecomposer decomposer;
        std::vector<int> node_ids(g.getNumNodes());
        std::iota(node_ids.begin(), node_ids.end(), 0);
        auto start = std::chrono::high_resolution_clock::now();
        TreeNode* DecompTree = decomposer.decompose(g, delta, node_ids);
        pure_oracle_running_times.push_back((std::chrono::duration<double, std::milli>( std::chrono::high_resolution_clock::now() - start)).count());

        computeNewDistances(g);
        g.print();
        return DecompTree;
    }


    void computeRLoads(int treeIndex, TreeNode*& tree, Graph& copyGraph) override {
            // --- 1️⃣ Ensure enough space for r-load maps ---
        if (m_idTree2edge2rload.size() <= treeIndex)
            m_idTree2edge2rload.resize(treeIndex + 1);

        auto& edge2Load = m_idTree2edge2rload[treeIndex];
        edge2Load.clear();

        // --- 2️⃣ BFS or DFS traversal (we'll use queue for clarity) ---
        std::queue<TreeNode*> q;
        q.push(tree);

        while (!q.empty()) {
            TreeNode* node = q.front();
            q.pop();

            // --- 3️⃣ Process each child: represents a cut S_child | V\S_child ---
            for (TreeNode* child : node->children) {
                // Add child to traversal queue
                q.push(child);

                const std::vector<int>& clusterVertices = child->members;
                if (clusterVertices.empty()) continue;

                // Build set for fast lookup
                std::unordered_set<int> S(clusterVertices.begin(), clusterVertices.end());

                // --- 4️⃣ Compute total cut capacity of this child cluster ---
                double cut = 0.0;
                for (int u : clusterVertices) {
                    for (int v : copyGraph.neighbors(u)) {
                        if (!S.count(v))  // boundary edge
                            cut += copyGraph.getEdgeCapacity(u, v);
                    }
                }
                if (cut <= 1e-12) cut = 1e-12;  // avoid zero-division

                // --- 5️⃣ Choose representative vertices (simple heuristic) ---
                int repParent = node->members.empty() ? clusterVertices[0] : node->members[0];
                int repChild  = clusterVertices[0];

                // --- 6️⃣ Compute shortest path between parent and child reps ---
                auto path = copyGraph.getShortestPath(repParent, repChild);
                if (path.size() < 2) continue;

                // --- 7️⃣ Update edge r-loads along that path ---
                for (size_t i = 0; i + 1 < path.size(); ++i) {
                    int u = path[i];
                    int v = path[i + 1];
                    double cap = copyGraph.getEdgeCapacity(u, v);
                    if (cap <= 1e-12) cap = 1e-12;

                    double delta = cut / cap;
                    edge2Load[{u, v}] += delta;
                    edge2Load[{v, u}] += delta;  // mirror for undirected graphs
                }
            }
        }
    }
};


using RaeckeCKRSolver = RaeckeFramework<CKRSolver, RaeckeCKRTransform>;


#endif //OBLIVIOUSROUTING_RAECKE_CKR_SOLVER_H