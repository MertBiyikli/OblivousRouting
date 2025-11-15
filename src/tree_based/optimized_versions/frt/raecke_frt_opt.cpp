//
// Created by Mert Biyikli on 05.11.25.
//
#include "raecke_frt_opt.h"

FRT_Tree RaeckeFRTOptimized::getTree(Graph_csr& g) {
    mcct.setGraph(g);
    mcct.setSeed(1337 + (int)m_trees.size()); // optional: vary per iteration
    mcct.setMaxPairsForUnitDemands(-1);       // all pairs, or cap if you want
    // If you have a custom demand set, add here via mcct.addDemand(s,t,w);

    FRT_Tree T = mcct.build(/*debug=*/false);

    computeNewDistance(g);
    return T;
}




void RaeckeFRTOptimized::setRequirements(const Graph_csr &g) {
    for (int u = 0; u < g.getNumNodes(); ++u) {
        for (int e = g.head[u]; e < g.head[u + 1]; ++e) {
            int v = g.edges[e].second;
            double cap = g.getEdgeCapacity(e);          // capacity, not distance
            if (cap <= 0) throw std::invalid_argument("Edge capacity must be positive.");
            mcct.addDemand(u, v, cap);                  // same as naïve: directed arc (u,v)
        }
    }
}




/*
 * TODO: Here lies sooo much potential for optimization.
 * 1) We can precompute the cut values for all nodes in the tree in a single pass.
 * 2) We can store parent pointers or use a more efficient data structure to
 *    avoid repeated shortest path computations.
 *
 */
void RaeckeFRTOptimized::computeRLoads(int treeIndex, FRT_Tree &_t, Graph_csr &copyGraph) {
    std::queue<std::shared_ptr<FRT_Node>> Q;
    for(auto& node : _t.GetRoot()->getChildren()) {
        Q.push(node);
    }

    if (m_tree2edge2rloads.size() <= treeIndex) {
        m_tree2edge2rloads.resize(treeIndex + 1);
    }

    m_tree2edge2rloads[treeIndex] = std::vector<double>(copyGraph.getNumEdges(), 0.0);
    auto& edge2Load = m_tree2edge2rloads[treeIndex];

    while (!Q.empty()) {
        auto node = Q.front();
        Q.pop();
        auto parent = node->getParent();
        if (parent->getCenter() == node->getCenter()) {
            for (auto& child : node->getChildren()) {
                Q.push(child);
            }
            continue;
        }
        std::vector<int> nodeVertices;
        for (const int& v: node->GetVertices()) {
            nodeVertices.push_back(v);
        }
        std::vector<int> remaining(copyGraph.getNumNodes());
        // fill with all vertices
        std::iota(remaining.begin(), remaining.end(), 0);
        // remove nodeVertices from remaining
        for (const int& v : nodeVertices) {
            auto it = std::find(remaining.begin(), remaining.end(), v);
            if (it != remaining.end()) {
                remaining.erase(it);
            }
        }

        double cut = 0.0;
        for(int u : nodeVertices) {
            for(const auto& [from, v] : copyGraph.neighbors(u)) {

                if(std::find(remaining.begin(), remaining.end(), v) != remaining.end()) {
                    cut += copyGraph.getEdgeCapacity(u, v);
                }
            }
        }

        int centerParent = parent->getCenter();
        int currentCenter = node->getCenter();

        auto path = copyGraph.getShortestPath(centerParent, currentCenter);
        for(int i = 0; i<path.size()-1; i++) {
            int u = path[i];
            int v = path[i+1];
            // if (u > v) std::swap(u, v); // Ensure (u,v) is in consistent order

            int eid = copyGraph.getEdgeId(u, v);
            int eid_rev = copyGraph.getEdgeId(v, u);

            if (eid < 0 || eid_rev < 0)
                throw std::runtime_error("Edge not found");

            edge2Load[eid]     += cut / copyGraph.getEdgeCapacity(eid);
            edge2Load[eid_rev] += cut / copyGraph.getEdgeCapacity(eid_rev);
        }
        for (auto& child : node->getChildren()) {
            Q.push(child);
        }
    }

    flow_store.begin(copyGraph.getNumEdges());  // ✅ resets buckets each iteration
    const auto demands = getAllDemands(500); // cap to 500 pairs if needed
    std::cout << "Demands count = " << demands.size() << "\n";

    std::vector<std::tuple<int,int,int,double>> local_flows;
    local_flows.reserve(demands.size() * 4);

    for (auto [s, t, d_val] : demands) {
        auto path_edges = copyGraph.getPathEdges(s, t);
        // print out edges
        for (int e : path_edges) {
            auto [u, v] = copyGraph.edges[e];
            std::cout << u << " -> " << v << " ";
        }
        for (int e : path_edges) {
            double flow_val = d_val;
            m_tree2edge2rloads[treeIndex][e] += flow_val;
            // print out added flow and its corresponding tree
            // std::cout << "\nAdding flow " << flow_val << " on edge " << e << " in tree " << treeIndex << "\n";
            local_flows.emplace_back(e, s, t, flow_val);
        }
    }

    for (auto [e, s, t, f] : local_flows)
        flow_store.add(e, s, t, f);


    flow_store.finalize();  // ✅ finalize only this tree’s flows
}
