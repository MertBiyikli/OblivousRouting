//
// Created by Mert Biyikli on 22.06.26.
//

#ifndef OBLIVIOUSROUTING_SOLVER_H
#define OBLIVIOUSROUTING_SOLVER_H

#include "routing/routing_table.h"
#include "data_structures/graph/Igraph.h"
#include "data_structures/graph/graph.h"
#include "errors.h"
#include <memory>

class ISolver {
protected:
    optimized::Graph<EdgeData>& graph;
public:
    explicit ISolver(optimized::Graph<EdgeData>& g) : graph(g) {}
    virtual ~ISolver() = default;

    virtual Result<std::unique_ptr<RoutingScheme>> solve() = 0;
    
    void resetEdgeDistance() {
        for (int v = 0; v < graph.getNumNodes(); ++v) {
            for (const auto& e : graph.edgesOf(v)) {
                auto& data = graph.edgeData(e.id);
                data.weight = 1.0 / graph.getNumUndirectedEdges();
            }
        }
    }
};
#endif //OBLIVIOUSROUTING_SOLVER_H