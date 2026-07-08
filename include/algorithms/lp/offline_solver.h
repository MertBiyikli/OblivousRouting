//
// Created by Mert Biyikli on 22.06.26.
//

#ifndef OBLIVIOUSROUTING_OFFLINE_SOLVER_H
#define OBLIVIOUSROUTING_OFFLINE_SOLVER_H

#include "core/solver.h"
#include "core/errors.h"
#include "utils/demands.h"
#include "routing/storage/allpair_routing_table.h"

class IOfflineSolver : public ISolver {

protected:
    demands demands_;
public:
    IOfflineSolver(IGraph& graph)
        : ISolver(graph) {}

    Result<std::unique_ptr<RoutingScheme>> solve() override {
        if (!demands_.size()) {
            throw std::runtime_error("[OfflineSolver]: demands must be initialized");
        }

        AllPairRoutingTable table;
        table.init(graph);
        auto flow = computeBasisFlows(table);
        if (!flow) {
            return getError(flow);
        }

        graph.resetEdgeDistance();

        return std::make_unique<AllPairRoutingScheme>(
            graph,
            std::move(table)
        );
    }

    virtual Result<void> computeBasisFlows(AllPairRoutingTable& table) = 0;

    void addDemand(int source, int target, double demand) {
        demands_.addDemand(source, target, demand);
    }

    double getDemandValue(int source, int target) const {
        auto value = demands_.getDemandValue(source, target);
        if (value.has_value()) {
            return value.value();
        } else {
            return 0.0;
        }
    }
};


#endif //OBLIVIOUSROUTING_OFFLINE_SOLVER_H