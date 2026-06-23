//
// Created by Mert Biyikli on 12.06.26.
//

#ifndef OBLIVIOUSROUTING_SEMI_ROUTING_ENGINE_H
#define OBLIVIOUSROUTING_SEMI_ROUTING_ENGINE_H

#include "../../data_structures/graph/Igraph.h"
#include "candidate_routing_scheme.h"

class IRoutingEngine {
public:
    virtual ~IRoutingEngine() = default;
    virtual const std::string getSolverBase() const = 0;

    virtual CandidateRoutingScheme preprocess(const IGraph& graph) = 0;
};

class ExistingSolverRoutingEngine final : public IRoutingEngine {
public:
    // TODO: in the future maybe make it also applicable for general oblivious routing as well
    explicit ExistingSolverRoutingEngine(std::shared_ptr<ILinearObliviousSolverBase> solver)
        : solver_(std::move(solver)) {}

    CandidateRoutingScheme preprocess(const IGraph& graph) override {
        CandidateRoutingScheme candidateScheme;

        // TODO:
        // 1. call existing solver
        // 2. obtain RoutingScheme
        // 3. decompose routing entries into candidate paths
        // 4. store them in CandidateRoutingScheme
        auto scheme = solver_->solve();
        for (int s = 0; s < graph.getNumNodes(); s++) {
            for (int t = 0; t<graph.getNumNodes(); t++) {
                if (s == solver_->getRootNode() || t == solver_->getRootNode()
                    || s == t) {
                    continue;
                }

                extractPath(graph, *scheme, s, t, candidateScheme);
            }
        }

        return candidateScheme;
    }

    void extractPath(const IGraph& g,const RoutingScheme& scheme,int s,int t,CandidateRoutingScheme& out) const {
        const int m = g.getNumDirectedEdges();

        std::vector<double> residual(m, 0.0);

        for (int e = 0; e < m; ++e) {
            const double f = scheme.getFlow(e, s, t);

            if (f > EPS) {
                residual[e] += f;
            } else if (f < -EPS) {
                const int anti = g.getAntiEdge(e);
                residual[anti] += -f;
            }
        }

        while (true) {
            std::vector<int> pathEdges;
            std::vector<bool> visited(g.getNumNodes(), false);

            const bool found = dfsDecompose(
                g,
                s,
                t,
                residual,
                pathEdges,
                visited
            );

            if (!found) break;

            double bottleneck = std::numeric_limits<double>::infinity();

            for (int e : pathEdges) {
                bottleneck = std::min(bottleneck, residual[e]);
            }

            if (bottleneck <= EPS || !std::isfinite(bottleneck)) {
                break;
            }

            for (int e : pathEdges) {
                residual[e] -= bottleneck;
                if (std::abs(residual[e]) < EPS) {
                    residual[e] = 0.0;
                }
            }

            Path path;
            path.source = s;
            path.target = t;
            path.edges = pathEdges;
            path.weight = bottleneck;

            out.addPath(s, t, std::move(path));
        }
    }

    virtual const std::string getSolverBase() const override{
        if (auto sol = dynamic_cast<ElectricalMWU*>(solver_.get())) {
            return "Electrical";
        } else if (auto sol = dynamic_cast<TreeMWU<FlatHST>*>(solver_.get())) {
            return "Tree";
        } else {
            return "UnknownSolverBase";
        }
    }

    bool dfsDecompose(const IGraph& g,int u,int t,std::vector<double>& residual,std::vector<int>& currentPath,std::vector<bool>& visited) const {
        if (u == t) {
            return true;
        }

        visited[u] = true;

        for (int neig: g.neighbors(u)) { // adapt name if needed
            int e = g.getEdgeId(u, neig);
            if (e == INVALID_EDGE_ID) continue;
            if (residual[e] <= EPS) continue;

            const auto [a, b] = g.getEdgeEndpoints(e);

            if (a != u) continue;

            const int v = b;

            if (visited[v]) continue;

            currentPath.push_back(e);

            if (dfsDecompose(g, v, t, residual, currentPath, visited)) {
                return true;
            }

            currentPath.pop_back();
        }

        return false;
    }

private:
    std::shared_ptr<ILinearObliviousSolverBase> solver_;
};

#endif //OBLIVIOUSROUTING_SEMI_ROUTING_ENGINE_H