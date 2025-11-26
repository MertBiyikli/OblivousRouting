//
// Created by Mert Biyikli on 23.10.25.
//

#include "ckr_partition.h"
#include <unordered_map>
#include <unordered_set>
#include "../../utils/priority_queue.h"
#include <random>

static constexpr int UNASSIGNED = -1;

void CKRPartition::init(const Graph &g, bool debug) {
    this->m_graph = g;
}

std::vector<int> CKRPartition::computePartition(const std::vector<int> &X_input, const double& Delta){
   const size_t n = m_graph.getNumNodes();
    const bool sparse_mode = (X_input.size() < 0.2 * n);  // heuristic switch

    // --- 1. Build random permutation π over X ---
    std::vector<int> X = X_input;
    std::mt19937_64 gen(std::random_device{}());
    std::shuffle(X.begin(), X.end(), gen);

    // --- 2. Sample R in [Δ/4, Δ/2] ---
    std::uniform_real_distribution<double> dist(Delta / 4.0, Delta / 2.0);
    const double R = dist(gen);

    // --- 3. Data structures ---
    // δ and P: store only what’s needed
    std::unordered_map<int, double> delta_map;       // for sparse mode
    std::unordered_map<int, int> label_map;

    std::vector<double> delta_dense;                 // for dense mode
    std::vector<int> P_dense;

    if (!sparse_mode) {
        delta_dense.assign(n, std::numeric_limits<double>::infinity());
        P_dense.assign(n, 0);
    } else {
        delta_map.reserve(X.size() * 4);
        label_map.reserve(X.size() * 4);
    }

    // Priority queue sized to what we’ll actually insert
    MinHeap<double, int> Q;
    Q.resize(X.size() * 4); // roughly |X|, not full n

    // --- 4. CKR loop ---
    for (int i = 0; i < (int)X.size(); ++i) {
        const int s = X[i];
        double& d_s = sparse_mode ? delta_map[s] : delta_dense[s];
        if (d_s > 0.0) {
            d_s = 0.0;
            Q.insertOrAdjustKey(s, 0.0);
        }

        while (!Q.empty() && Q.topKey() <= R) {
            const int w = Q.top();
            const double dw = Q.topKey();
            Q.deleteTop();

            // label assignment
            if (sparse_mode) {
                if (!label_map.count(w)) label_map[w] = i + 1;
            } else {
                if (P_dense[w] == 0) P_dense[w] = i + 1;
            }

            // relax neighbors
            for (int u : m_graph.neighbors(w)) {
                const double nd = dw + m_graph.getEdgeDistance(w, u);
                double& du = sparse_mode ? delta_map[u] : delta_dense[u];
                if (nd < du) {
                    du = nd;
                    Q.insertOrAdjustKey(u, nd);
                }
            }
        }
    }

    // --- 5. Collect labels for X only ---
    std::vector<int> labels(X.size(), UNASSIGNED);
    for (size_t k = 0; k < X.size(); ++k) {
        const int v = X[k];
        labels[k] = sparse_mode
                        ? (label_map.count(v) ? label_map[v] - 1 : UNASSIGNED)
                        : (P_dense[v] ? P_dense[v] - 1 : UNASSIGNED);
    }
    return labels;
}


 // namespace MendelScaling
std::vector<int> CKRPartition::computePartition(const std::vector<int>& X_input, const double& Delta, CKRLevel& L) {
    const size_t n = m_graph.getNumNodes();
    const bool sparse_mode = (X_input.size() < 0.2 * n);  // heuristic switch

    // --- 1. Build random permutation π over X ---
    std::vector<int> X = X_input;
    std::mt19937_64 gen(std::random_device{}());
    std::shuffle(X.begin(), X.end(), gen);

    // --- 2. Sample R in [Δ/4, Δ/2] ---
    std::uniform_real_distribution<double> dist(Delta / 4.0, Delta / 2.0);
    const double R = dist(gen);

    L.R = R;
    L.owner.assign(n, -1);
    L.pred.assign(n, -1);

    L.centers.clear();
    L.centers.reserve(n);

    // --- 3. Data structures ---
    // δ and P: store only what’s needed
    /*
    std::unordered_map<int, double> delta_map;       // for sparse mode
    std::unordered_map<int, int> label_map;
*/
    std::vector<double> delta_dense;                 // for dense mode
    std::vector<int> P_dense;

    delta_dense.assign(n, std::numeric_limits<double>::infinity());
    P_dense.assign(n, 0);
    /*
    if (!sparse_mode) {
        delta_dense.assign(n, std::numeric_limits<double>::infinity());
        P_dense.assign(n, 0);
    } else {
        delta_map.reserve(X.size() * 4);
        label_map.reserve(X.size() * 4);
    }*/




    // Priority queue sized to what we’ll actually insert
    MinHeap<double, int> Q;
    Q.resize(X.size() * 4); // roughly |X|, not full n

    // --- 4. CKR loop ---
    for (int i = 0; i < (int)X.size(); ++i) {
        const int s = X[i];

        // if node has been assigned already, skip
        if (L.owner[s]!=-1) continue; // already assigned
        // otherwise make it to a new center
        L.centers.push_back(s); // add center
        L.pred[s] = -1; // center has no predecessor

        double& d_s = delta_dense[s];
        if (d_s > 0.0) {
            d_s = 0.0;
            Q.insertOrAdjustKey(s, 0.0);
        }

        while (!Q.empty() && Q.topKey() <= R) {
            const int w = Q.top();
            const double dw = Q.topKey();
            Q.deleteTop();

            if (L.owner[w]!=-1) continue; // already assigned

            L.owner[w] = s; // assign owner

            // label assignment
            if (P_dense[w] == 0) P_dense[w] = i + 1;

            // relax neighbors
            for (int u : m_graph.neighbors(w)) {
                if (L.owner[u]!=-1) continue; // already assigned
                const double nd = dw + m_graph.getEdgeDistance(w, u);
                double& du = delta_dense[u];
                // double& du = sparse_mode ? delta_map[u] : delta_dense[u];
                if (nd < du) {
                    du = nd;
                    Q.insertOrAdjustKey(u, nd);
                    L.pred[u] = w; // set predecessor
                }
            }
        }
    }

    // --- 5. Collect labels for X only ---
    std::vector<int> labels(X.size(), UNASSIGNED);
    for (size_t k = 0; k < X.size(); ++k) {
        const int v = X[k];
        labels[k] = (P_dense[v] ? P_dense[v] - 1 : UNASSIGNED);
        /*
        labels[k] = sparse_mode
                        ? (label_map.count(v) ? label_map[v] - 1 : UNASSIGNED)
                        : (P_dense[v] ? P_dense[v] - 1 : UNASSIGNED);
        */
    }
    return labels;
}



 // namespace MendelScaling
std::vector<int> CKRPartition::computePartition(const Graph_csr& g, const std::vector<int>& X_input, const double& Delta, CKRLevel& L) {
    const size_t n = g.getNumNodes();


    // --- 1. Build random permutation π over X ---
    std::vector<int> X = X_input;
    std::mt19937_64 gen(std::random_device{}());
    std::shuffle(X.begin(), X.end(), gen);

    // --- 2. Sample R in [Δ/4, Δ/2] ---
    std::uniform_real_distribution<double> dist(Delta / 4.0, Delta / 2.0);
    const double R = dist(gen);

    L.R = R;
    L.owner.assign(n, -1);
    L.pred.assign(n, -1);

    L.centers.clear();
    L.centers.reserve(n);

    // --- 3. Data structures ---
    // δ and P: store only what’s needed
    std::vector<double> delta_dense;                 // for dense mode
    std::vector<int> P_dense;

    delta_dense.assign(n, std::numeric_limits<double>::infinity());
    P_dense.assign(n, 0);




    // Priority queue sized to what we’ll actually insert
    MinHeap<double, int> Q;
    Q.resize(X.size() * 4); // roughly |X|, not full n

    // --- 4. CKR loop ---
    for (int i = 0; i < (int)X.size(); ++i) {
        const int s = X[i];

        // if node has been assigned already, skip
        if (L.owner[s]!=-1) continue; // already assigned
        // otherwise make it to a new center
        L.centers.push_back(s); // add center
        L.pred[s] = -1; // center has no predecessor

        double& d_s = delta_dense[s];
        if (d_s > 0.0) {
            d_s = 0.0;
            Q.insertOrAdjustKey(s, 0.0);
        }

        while (!Q.empty() && Q.topKey() <= R) {
            const int w = Q.top();
            const double dw = Q.topKey();
            Q.deleteTop();

            if (L.owner[w]!=-1) continue; // already assigned

            L.owner[w] = s; // assign owner

            // label assignment
            if (P_dense[w] == 0) P_dense[w] = i + 1;

            // relax neighbors
            for (auto& u : g.neighbors(w)) {
                if (L.owner[u]!=-1) continue; // already assigned
                const double nd = dw + g.getEdgeDistance(w, u);
                double& du = delta_dense[u];
                // double& du = sparse_mode ? delta_map[u] : delta_dense[u];
                if (nd < du) {
                    du = nd;
                    Q.insertOrAdjustKey(u, nd);
                    L.pred[u] = w; // set predecessor
                }
            }
        }
    }

    // --- 5. Collect labels for X only ---
    std::vector<int> labels(X.size(), UNASSIGNED);
    for (size_t k = 0; k < X.size(); ++k) {
        const int v = X[k];
        labels[k] = (P_dense[v] ? P_dense[v] - 1 : UNASSIGNED);

    }
    return labels;
}
