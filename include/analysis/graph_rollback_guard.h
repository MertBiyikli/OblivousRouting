//
// Created by Mert on 28.10.26.
//

#ifndef E_ROUTING_GRAPH_ROLLBACK_GUARD_H
#define E_ROUTING_GRAPH_ROLLBACK_GUARD_H

#include "data_structures/graph/graph.h"

// TODO: move this potential to graph.h
class GraphRollbackGuard {
  public:
    explicit GraphRollbackGuard(optimized::Graph<EdgeData>& graph) : graph_(graph), checkpoint_(graph.checkpoint()) {}

    ~GraphRollbackGuard() noexcept {
        if (active_) {
            auto result = graph_.rollback(checkpoint_);

            if (!result) {
                std::terminate();
            }
        }
    }

    GraphRollbackGuard(const GraphRollbackGuard&) = delete;
    GraphRollbackGuard& operator=(const GraphRollbackGuard&) = delete;

    int checkpoint() const noexcept { return checkpoint_; }

  private:
    optimized::Graph<EdgeData>& graph_;
    int checkpoint_;
    bool active_ = true;
};

#endif // E_ROUTING_GRAPH_ROLLBACK_GUARD_H
