//
// Created by Mert Biyikli on 29.11.25.
//

#ifndef OBLIVIOUSROUTING_EFFICIENT_RAECKE_BASE_H
#define OBLIVIOUSROUTING_EFFICIENT_RAECKE_BASE_H


#include "../datastructures/IGraph.h"
#include "raecke_oracle_iteration.h"
#include <cassert>
#include <chrono>
#include <cmath>

template<typename Tree, typename DistanceType>
class EfficientRaeckeBase {
public:
    double lambdaSum = 0;
    double diameter = 0;
    IGraph* g_ptr = nullptr;

    std::vector<double> rload_current;
    std::vector<double> rload_total;

    std::vector<OracleTreeIteration> iteration;
    std::vector<double> oracle_running_times;
    std::vector<double> pure_oracle_running_times;

    virtual Tree getTree() = 0;


    void run() {
        if (!g_ptr)
            throw std::runtime_error("CKR: graph pointer is null");

        lambdaSum = 0.0;
        while (lambdaSum < 1.0) {
            auto start = std::chrono::high_resolution_clock::now();
            lambdaSum += iterate();

            oracle_running_times.push_back(
                std::chrono::duration<double, std::milli>(
                    std::chrono::high_resolution_clock::now() - start
                ).count()
            );

        }

        g_ptr = nullptr;
    }


    double iterate() {
        Tree t = getTree();
        assert(t != nullptr);
        computeRLoads(t);
        double l = getMaxRload();
        double lambda = std::min(1.0/l, 1.0 - lambdaSum);

        // update weights
        lambdaSum += lambda;
        computeNewDistances(lambda);

        OracleTreeIteration iter(lambda, t, g_ptr->getDistanceVector<DistanceType>());

        iteration.push_back(iter);

        return lambda;
    }
    void setGraph(IGraph& graph) {
        g_ptr = &graph;              // no copy
        init();
    }

    void init() {
        if (!g_ptr) {
            throw std::runtime_error("No given CKR graph");
        }
        g_ptr->finalize();
        diameter = g_ptr->GetDiameter();

        rload_current.assign(g_ptr->getNumEdges(), 0.0);
        rload_total.assign(g_ptr->getNumEdges(), 0.0);
    }


    void addCurrentLoad(double lambda) {
        assert(g_ptr != nullptr);
        for (int e = 0; e < g_ptr->getNumEdges(); ++e) {
            rload_total[e] += rload_current[e] * lambda;
        }
    }

    void computeRLoads(Tree t) {
        assert(t != nullptr);
        assert(g_ptr != nullptr);
        // clear current rloads

        std::fill(rload_current.begin(), rload_current.end(), 0.0);


        std::queue<Tree> q;
        q.push(t);
        while (!q.empty()) {
            Tree node = q.front();
            q.pop();

            // --- 3️⃣ Process each child: represents a cut S_child | V\S_child ---
            for (auto child : node->children) {
                // Add child to traversal queue
                q.push(child);


                // if node has the same nodes as child, skip
                if (node->members.size() == child->members.size()) {
                    continue;
                }
                const std::vector<int>& clusterVertices = child->members;
                if (clusterVertices.empty()) continue;

                // Build set for fast lookup
                std::vector<char> S(g_ptr->getNumNodes(), 0);
                for (int v : clusterVertices) S[v] = 1;


                // --- 4️⃣ Compute total cut capacity of this child cluster ---
                double cut = 0.0;
                for (int u : clusterVertices) {
                    for (auto&  v : g_ptr->neighbors(u)) {
                        if (!S[v])  // boundary edge
                            cut += g_ptr->getEdgeCapacity(u, v);
                    }
                }
                if (cut <= 1e-12) cut = 1e-12;  // avoid zero-division

                // --- 5️⃣ Choose representative vertices (simple heuristic) ---
                int repParent = node->members.empty() ? clusterVertices[0] : node->members[0];
                int repChild  = clusterVertices[0];

                // --- 6️⃣ Compute shortest path between parent and child reps ---
                auto path = g_ptr->getShortestPath(repParent, repChild);
                if (path.size() < 2) continue;

                // --- 7️⃣ Update edge r-loads along that path ---
                for (size_t i = 0; i + 1 < path.size(); ++i) {
                    int u = path[i];
                    int v = path[i + 1];

                    int e = g_ptr->getEdgeId(u, v);
                    double rload = rload_current[e];
                    double cap = g_ptr->getEdgeCapacity(u, v);
                    if (cap <= 1e-12) cap = 1e-12;

                    rload += cut / cap;
                    rload_current[e] = rload;

                    int anti_e = g_ptr->getEdgeId(v, u);
                    rload_current[anti_e] = rload;

                }
            }
        }
    }
    double getMaxRload() const {
        assert(rload_current.size() > 0);

        double max_load = 0.0;
        for (const auto& load : rload_current) {
            if (load > max_load) {
                max_load = load;
            }
        }
        return max_load;
    }

    void computeNewDistances(double lambda) {
        assert(g_ptr != nullptr);
        constexpr double EPS = 1e-12;

        addCurrentLoad(lambda);

        double max_r = 0.0;
        for (auto& r : rload_total) max_r = std::max(max_r, r);

        double sumExp = 0.0;
        for (auto& r : rload_total) sumExp += std::exp(r - max_r);
        if (sumExp <= 0.0 || !std::isfinite(sumExp)) sumExp = 1.0;

        // 3) distances
        double min_d = std::numeric_limits<double>::infinity();
        std::vector<double> newDist;
        newDist.reserve(rload_total.size());

        for (int e = 0; e < g_ptr->getNumEdges(); ++e) {
            double r = rload_total[e];

            double cap = g_ptr->getEdgeCapacity(e);
            if (cap < EPS) cap = EPS;
            double d = (std::exp(r - max_r) / cap) / sumExp;
            newDist[e] = d;
            if (d < min_d) min_d = d;
        }
        if (min_d < EPS) min_d = EPS;

        for (int e = 0; e < g_ptr->getNumEdges(); ++e) {
            double d = newDist[e];
            double norm = d / min_d;
            if (norm < 1.0) norm = 1.0;
            g_ptr->updateEdgeDistance(e, norm);
        }
    }

    const int getIterationCount() const {
        return static_cast<int>(oracle_running_times.size());
    }


    const std::vector<double> getOracleRunningTimes() const {
        return oracle_running_times;
    }

    const std::vector<double> getPureOracleRunningTimes() const {
        return pure_oracle_running_times;
    }


    std::vector<OracleTreeIteration>& getIterations() {
        return iteration;
    }

};

#endif //OBLIVIOUSROUTING_EFFICIENT_RAECKE_BASE_H