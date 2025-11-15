//
// Created by Mert Biyikli on 27.10.25.
//

#include "raecke_ckr_optimized.h"

#include <cassert>

void RaeckeCKR::init(Graph &g) {
    this->m_graph = g;
    for (int i = 0; i < g.getNumNodes(); ++i) {
        for (auto& u : g.neighbors(i)) {
            if ( i < u) {
                double cap = g.getEdgeCapacity(i, u);
                m_graph.updateEdgeDistance(i, u, cap);
            }
        }
    }
}

TreeNode* RaeckeCKR::getTree(Graph& g)
{
    const int n = g.getNumNodes();
    if (n == 0) return nullptr;

    // ----- 1) choose base scale -----
    if (!g.IsDistanceMatrixComputed())
        g.createDistanceMatrix();

    double Delta = g.GetDiameter();

    std::mt19937 rng(std::random_device{}());
    m_levels.clear();
    m_levels.reserve(32);

    // Finest level: each vertex is a center
    std::vector<int> current_centers(n);
    std::iota(current_centers.begin(), current_centers.end(), 0);

    // For assembling the HST tree
    std::vector<TreeNode*> node_ptrs(n, nullptr);
    for (int v = 0; v < n; ++v) {
        node_ptrs[v] = new TreeNode(v); // leaf nodes
        node_ptrs[v]->members.push_back(v);
    }

    TreeNode* root = nullptr;

    // ----- 2) hierarchical CKR levels -----
    while (current_centers.size() > 1) {
        m_levels.emplace_back();
        CKRLevel& L = m_levels.back();
        build_ckr_level(g, Delta, L, rng);

        // ----- build clusters -> TreeNodes -----
        std::unordered_map<int, TreeNode*> center_to_node;
        center_to_node.reserve(L.centers.size()*2);

        // group vertices by owner
        for (int v = 0; v < n; ++v) {
            int c = L.owner[v];
            if (c == -1) continue;
            TreeNode* parent;
            auto it = center_to_node.find(c);
            if (it == center_to_node.end()) {
                parent = new TreeNode(c);
                parent->radius = Delta;
                center_to_node[c] = parent;
            } else parent = it->second;

            parent->children.push_back(node_ptrs[v]);
            node_ptrs[v]->parent = parent;
            // merge clusters (this is key!)
            parent->members.insert(parent->members.end(),
                                   node_ptrs[v]->members.begin(),
                                   node_ptrs[v]->members.end());
        }

        // prepare for next level
        node_ptrs.clear();
        node_ptrs.reserve(center_to_node.size());
        for (auto& [c,node] : center_to_node)
            node_ptrs.push_back(node);

        current_centers = L.centers;
        Delta /= 2;  // shrink scale for next level
    }

    // ----- 3) final root -----
    if (!node_ptrs.empty()) root = node_ptrs.front();
    if (root) root->parent = nullptr;

    assert(root && (int)root->members.size() == n);
    return root;
}


/*
TreeNode* RaeckeCKR::getTree(Graph &g) {
    m_graph.createDistanceMatrix();
    double delta = m_graph.GetDiameter();
    TreeDecomposer decomposer;
    std::vector<int> node_ids(g.getNumNodes());
    std::iota(node_ids.begin(), node_ids.end(), 0);
    return decomposer.decompose(g, delta, node_ids);
}
*/


void RaeckeCKR::computeRLoads(int idx, TreeNode *t, Graph &g) {
    // iterate through the tree and store the rloads into the adjacency list
    std::queue<TreeNode*> q;
        q.push(t);

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
                std::vector<char> S(m_graph.getNumNodes(), 0);
                for (int v : clusterVertices) S[v] = 1;


                // --- 4️⃣ Compute total cut capacity of this child cluster ---
                double cut = 0.0;
                for (int u : clusterVertices) {
                    for (int v : m_graph.neighbors(u)) {
                        if (!S[v])  // boundary edge
                            cut += m_graph.getEdgeCapacity(u, v);
                    }
                }
                if (cut <= 1e-12) cut = 1e-12;  // avoid zero-division

                // --- 5️⃣ Choose representative vertices (simple heuristic) ---
                int repParent = node->members.empty() ? clusterVertices[0] : node->members[0];
                int repChild  = clusterVertices[0];

                // --- 6️⃣ Compute shortest path between parent and child reps ---
                auto path = m_graph.getShortestPath(repParent, repChild);
                if (path.size() < 2) continue;

                // --- 7️⃣ Update edge r-loads along that path ---
                for (size_t i = 0; i + 1 < path.size(); ++i) {
                    int u = path[i];
                    int v = path[i + 1];
                    double cap = m_graph.getEdgeCapacity(u, v);
                    if (cap <= 1e-12) cap = 1e-12;

                    double delta = cut / cap;
                    addLoadToEdge(u, v, delta);
                    /*
                    edge2Load[{u, v}] += delta;
                    edge2Load[{v, u}] += delta;  // mirror for undirected graphs*/
                }
            }
        }
}


// TODO: this is for now a naive not optimized, for each newly added load => O(n) linear possible -> could be reduced to O(log) at least with binary search tree
void RaeckeCKR::addLoadToEdge(int u, int v, double load) {
    edge2Load[{u, v}] += load;
    edge2Load[{v, u}] += load;
    /*
    bool found_forward = false, found_backwards = false;
    // find the edge in the adjacency list and add the load
    if (m_adj_edge_loads.size() <= u) {
        m_adj_edge_loads.resize(u + 1);
        m_edgeLoads.resize(u + 1);
    }

    if (m_adj_edge_loads.size() <= v) {
        m_adj_edge_loads.resize(v + 1);
        m_edgeLoads.resize(v + 1);
    }

    for (size_t i = 0; i < m_adj_edge_loads[u].size(); i++) {
        if (m_adj_edge_loads[u][i] == v) {
            m_edgeLoads[u][i] += load;
            found_forward = true;
        }
    }

    // also for the reverse edge
    for (size_t i = 0; i < m_adj_edge_loads[v].size(); i++) {
        if (m_adj_edge_loads[v][i] == u) {
            m_edgeLoads[v][i] += load;
            found_backwards = true;
        }
    }

    if ( !found_forward ) {
        // add flow to the adjacency list
        m_adj_edge_loads[u].push_back(v);
        m_edgeLoads[u].push_back(load);
    }

    if (!found_backwards) {
        // add flow to the adjacency list
        m_adj_edge_loads[v].push_back(u);
        m_edgeLoads[v].push_back(load);
    }*/
}


double RaeckeCKR::getMaxRload() {
    double max_load = 0.0;

    for (const auto& [edge, load] : edge2Load) {
        if (load > max_load) {
            max_load = load;
        }
    }
    /*
    for (int u = 0; u < m_graph.getNumNodes(); ++u) {
        for (size_t i = 0; i < m_adj_edge_loads[u].size(); ++i) {
            double load = m_edgeLoads[u][i];
            if (load > max_load) {
                max_load = load;
            }
        }
    }*/
    return max_load;
}

void RaeckeCKR::computeNewDistances(Graph &g) {
    constexpr double EPS = 1e-12;

    // 1) Build total_r(e) for all edges, default 0
    std::unordered_map<std::pair<int,int>, double> total_r;
    total_r.reserve(g.getNumEdges());

    for (int u = 0; u < g.getNumNodes(); ++u) {
        for (int v : g.neighbors(u)) {
            int a = std::min(u,v), b = std::max(u,v);
            // only one direction
            if (u != a) continue;
            auto it = edge2Load.find({a,b});
            total_r[{a,b}] = (it == edge2Load.end()) ? 0.0 : it->second;
        }
    }

    // 2) softmax with shift
    double max_r = 0.0;
    for (auto& [e,r] : total_r) max_r = std::max(max_r, r);

    double sumExp = 0.0;
    for (auto& [e,r] : total_r) sumExp += std::exp(r - max_r);
    if (sumExp <= 0.0 || !std::isfinite(sumExp)) sumExp = 1.0;

    // 3) distances
    double min_d = std::numeric_limits<double>::infinity();
    std::unordered_map<std::pair<int,int>, double> newDist;
    newDist.reserve(total_r.size());

    for (auto& [e,r] : total_r) {
        auto [a,b] = e;
        double cap = g.getEdgeCapacity(a,b);
        if (cap < EPS) cap = EPS;
        double d = (std::exp(r - max_r) / cap) / sumExp;
        newDist[e] = d;
        if (d < min_d) min_d = d;
    }
    if (min_d < EPS) min_d = EPS;

    for (auto& [e,d] : newDist) {
        auto [a,b] = e;
        double norm = d / min_d;
        if (norm < 1.0) norm = 1.0;
        g.updateEdgeDistance(a,b, norm);
        g.updateEdgeDistance(b,a, norm); // mirror, if your graph stores both arcs
    }
}


