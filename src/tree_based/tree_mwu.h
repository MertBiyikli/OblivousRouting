//
// Created by Mert Biyikli on 09.02.26.
//

#ifndef OBLIVIOUSROUTING_TREE_MWU_H
#define OBLIVIOUSROUTING_TREE_MWU_H

#include "../solver/solver.h"
#include "../solver/mwu_framework.h"
#include "raecke_oracle_iteration.h"
#include "raecke_transform.h"
#include "hst.h"
#include "tree_oracle.h"
#include "tree_transform.h"


class TreeMWU : public AllPairObliviousSolverBase, public MWUFramework {

    std::unique_ptr<TreeOracle> oracle;

    // member for storing the information for the computation
    double lambda_sum = 0.0;
    std::vector<double> rload_current, rload_total;
    std::vector<double> current_distances;
    std::vector<TreeIteration> iteration;
public:
    TreeMWU(IGraph& g, std::unique_ptr<TreeOracle> _oracle) : AllPairObliviousSolverBase(g), oracle(std::move(_oracle)) {
        current_distances.assign(graph.getNumEdges(), 1.0);
        rload_current.assign(graph.getNumEdges(), 0.0);
        rload_total.assign(graph.getNumEdges(), 0.0);
    }


    virtual void computeBasisFlows(AllPairRoutingTable &table) override {
        auto t0 = timeNow();
        run();
        this->solve_time = duration((timeNow()-t0));
        normalizeLambdas();

        t0 = timeNow();
        transformSolution(table);
        this->transformation_time = duration((timeNow()-t0));
    }

    virtual void transformSolution(AllPairRoutingTable& table) {
        TreeTransform transform(graph, iteration);
        transform.transform();
        table = transform.getRoutingTable();

        //table.printFlows(graph);
    }


    // MWU loop
    void run() {
        lambda_sum = 0.0;
        while (lambda_sum < 1.0) {
            lambda_sum += treeOracle();
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
        double lambda = std::min(1.0/l, 1.0 - lambda_sum);

        TreeIteration iter( std::move(t), current_distances, lambda);
        iteration.push_back(iter);

        computeNewDistances(lambda);

        return lambda;
    }

    void normalizeLambdas() {
        for (auto& it : iteration) {
            double lambda = it.getLambda();
            lambda /= lambda_sum;
            it.setLambda(lambda);
        }
    }


    void addCurrentLoad(double lambda) {
        for (int e = 0; e < graph.getNumEdges(); ++e) {
            rload_total[e] += rload_current[e] * lambda;
        }
    }

    void computeRLoads(const std::shared_ptr<HSTNode> t) {
        assert(t != nullptr);
        // clear current rloads

        std::fill(rload_current.begin(), rload_current.end(), 0.0);



        std::queue<std::shared_ptr<HSTNode>> q;
        q.push(t);
        while (!q.empty()) {
            auto node = q.front();
            q.pop();

            // Process each child: represents a cut S_child | V\S_child ---
            for (const auto& child : node->getChildren()) {
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


                // Compute total cut capacity of this child cluster
                double cut = 0.0;
                for (int u : clusterVertices) {
                    for (auto&  v : graph.neighbors(u)) {
                        if (!S[v])  // boundary edge
                            cut += graph.getEdgeCapacity(u, v);
                    }
                }
                if (cut <= 1e-12) cut = 1e-12;  // avoid zero-division

                // Choose representative vertices
                int repParent = node->getMembers().empty() ? clusterVertices[0] : node->getMembers()[0];
                int repChild  = clusterVertices[0];

                // Compute shortest path between parent and child reps
                auto path = graph.getShortestPath(repParent, repChild);
                if (path.size() < 2) continue;

                // Update edge r-loads along that path
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

    void activateMendelScaling(bool flag) {
        if (oracle) {
            oracle->applyMendelScaling = flag;
        }else {
            throw std::runtime_error("TreeMWU: Oracle must be set before activating Mendel Scaling.");
        }
    }
};

#endif //OBLIVIOUSROUTING_TREE_MWU_H