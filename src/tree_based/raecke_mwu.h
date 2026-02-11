//
// Created by Mert Biyikli on 08.12.25.
//

#ifndef OBLIVIOUSROUTING_RAECKE_LINEAR_OR_H
#define OBLIVIOUSROUTING_RAECKE_LINEAR_OR_H


#include "raecke_oracle_iteration.h"
#include "../solver/solver.h"
#include "../solver/mwu_framework.h"
#include "raecke_oracle.h"
#include "raecke_oracle_iteration.h"
#include "raecke_tree.h"
#include "raecke_transform.h"

class RaeckeMWU : public AllPairObliviousSolverBase, public MWUFramework {
public:
    std::unique_ptr<RaeckeOracle> oracle;
    double lambdaSum;

    std::vector<double> rload_current;
    std::vector<double> rload_total;

    std::vector<OracleTreeIteration> iteration;

    std::vector<double> current_distances;

    RaeckeMWU(IGraph& g,
              std::unique_ptr<RaeckeOracle> _oracle):AllPairObliviousSolverBase(g),
    oracle(std::move(_oracle)),
          lambdaSum(0.0) {

        current_distances.assign(graph.getNumEdges(), 1.0);
        rload_current.assign(graph.getNumEdges(), 0.0);
        rload_total.assign(graph.getNumEdges(), 0.0);
    }


    virtual void transformSolution(AllPairRoutingTable& table) {
        EfficientRaeckeTransform transform(graph, iteration);
        transform.transform();
        table = transform.getRoutingTable();
        // table.printFlows(graph);

    }

    void computeBasisFlows(AllPairRoutingTable &table) override {
        auto t0 = timeNow();
        run();
        this->solve_time = duration((timeNow()-t0));
        normalizeLambdas();
        transformSolution(table);
        //table.printFlows(graph);

    }

    void run() {
        lambdaSum = 0.0;
        while (lambdaSum < 1.0) {
            lambdaSum += treeOracle();
            iteration_count++;
        }
    }

    double treeOracle() {
        auto start = timeNow();
        auto t = oracle->getTree(current_distances);
        this->oracle_running_times.push_back(duration((timeNow()-start)));

        assert(t != nullptr);
        computeRLoads(t);
        double l = getMaxRload();
        double lambda = std::min(1.0/l, 1.0 - lambdaSum);

        OracleTreeIteration iter(lambda, t, current_distances);
        iteration.push_back(iter);

        computeNewDistances(lambda);

        return lambda;
    }

    void normalizeLambdas() {
        for (auto& it : iteration) {
            double lambda = it.getLambda();
            lambda /= lambdaSum;
            it.setLambda(lambda);
        }
    }


    void addCurrentLoad(double lambda) {
        for (int e = 0; e < graph.getNumEdges(); ++e) {
            rload_total[e] += rload_current[e] * lambda;
        }
    }

    void computeRLoads(std::shared_ptr<ITreeNode> t) {
        assert(t != nullptr);
        // clear current rloads

        std::fill(rload_current.begin(), rload_current.end(), 0.0);



        std::queue<std::shared_ptr<ITreeNode>> q;
        q.push(t);
        while (!q.empty()) {
            auto node = q.front();
            q.pop();

            // --- 3️⃣ Process each child: represents a cut S_child | V\S_child ---
            for (auto child : node->getChildren()) {
                // Add child to traversal queue
                q.push(child);


                // if node has the same nodes as child, skip
                if (node->getMembers().size() == child->getMembers().size()) {
                    continue;
                }
                const std::vector<int>& clusterVertices = child->getMembers();
                if (clusterVertices.empty()) continue;

                // Build set for fast lookup
                std::vector<char> S(graph.getNumNodes(), 0);
                for (int v : clusterVertices) S[v] = 1;


                // --- 4️⃣ Compute total cut capacity of this child cluster ---
                double cut = 0.0;
                for (int u : clusterVertices) {
                    for (auto&  v : graph.neighbors(u)) {
                        if (!S[v])  // boundary edge
                            cut += graph.getEdgeCapacity(u, v);
                    }
                }
                if (cut <= 1e-12) cut = 1e-12;  // avoid zero-division

                // --- 5️⃣ Choose representative vertices (simple heuristic) ---
                int repParent = node->getMembers().empty() ? clusterVertices[0] : node->getMembers()[0];
                int repChild  = clusterVertices[0];

                // --- 6️⃣ Compute shortest path between parent and child reps ---
                auto path = graph.getShortestPath(repParent, repChild);
                if (path.size() < 2) continue;

                // --- 7️⃣ Update edge r-loads along that path ---
                for (size_t i = 0; i + 1 < path.size(); ++i) {
                    int u = path[i];
                    int v = path[i + 1];

                    int e = graph.getEdgeId(u, v);
                    double rload = rload_current[e];
                    double cap = graph.getEdgeCapacity(u, v);
                    if (cap <= 1e-12) cap = 1e-12;

                    rload += cut / cap;
                    rload_current[e] = rload;

                    int anti_e = graph.getEdgeId(v, u);
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

        for (int e = 0; e < graph.getNumEdges(); ++e) {
            double r = rload_total[e];

            double cap = graph.getEdgeCapacity(e);
            if (cap < EPS) cap = EPS;
            double d = (std::exp(r - max_r) / cap) / sumExp;
            newDist[e] = d;
            if (d < min_d) min_d = d;
        }
        if (min_d < EPS) min_d = EPS;

        for (int e = 0; e < graph.getNumEdges(); ++e) {
            double d = newDist[e];
            double norm = d / min_d;
            if (norm < 1.0) norm = 1.0;
            graph.updateEdgeDistance(e, norm);
            current_distances[e] = norm;
        }
    }

};



#endif //OBLIVIOUSROUTING_RAECKE_LINEAR_OR_H