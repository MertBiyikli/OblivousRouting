#pragma once

#include "mwu_framework.h"

class FlowSparsifier : public MWUFramework{
    public:
    FlowSparsifier(IGraph& g, int root):
    MWUFramework(g, root){};

    std::vector<double> _mwu_weights;
    std::vector<double> _rel_loads;
    std::vector<LinearRoutingTable> tables;


    virtual Result<void> computeBasisFlows(LinearRoutingTable& table) override {
        table.init(graph);
        auto res = run(table);
        if (!res) {
            return getError(res);
        }

        res = scaleFlowDown(table);
        if (!res) {
            return getError(res);
        }

        return {};
    }

    virtual Result<void> updateDistances(const std::vector<double> &distances) override {
        for (int e = 0; e < graph.getNumDirectedEdges(); ++e) {
            _mwu_weights[e] =  distances[e];
        }
        return {};
    }

    virtual Result<void> run(LinearRoutingTable &table);
    virtual Result<void> scaleFlowDown(LinearRoutingTable &table);
    virtual void printAdditionalStats() {
        return;
    }


    void computeNewDistances() {
        std::vector<double> current_distances(graph.getNumDirectedEdges());
        double max_r = 0.0;
        for (auto& r : _rel_loads) max_r = std::max(max_r, r);

        double sumExp = 0.0;
        for (auto& r : _rel_loads) sumExp += std::exp(r - max_r);
        if (sumExp <= 0.0 || !std::isfinite(sumExp)) sumExp = 1.0;

        double min_d = std::numeric_limits<double>::infinity();
        std::vector<double> newDist;
        newDist.resize(_rel_loads.size());

        for (int e = 0; e < graph.getNumDirectedEdges(); ++e) {
            double r = _rel_loads[e];

            double cap = graph.getEdgeCapacity(e);
            if (cap < EPS) cap = EPS;
            double d = (std::exp(r - max_r) / cap) / sumExp;
            newDist[e] = d;
            if (d < min_d) min_d = d;
        }
        if (min_d < EPS) min_d = EPS;

        for (int e = 0; e < graph.getNumDirectedEdges(); ++e) {
            double d = newDist[e];
            double norm = d / min_d;
            if (norm < 1.0) norm = 1.0;
            current_distances[e] = norm;
        }
        updateDistances(current_distances);
    }
};