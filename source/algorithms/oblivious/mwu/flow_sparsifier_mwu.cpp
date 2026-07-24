//
// Created by Mert Biyikli on 21.07.26.
//

#include "algorithms/oblivious/mwu/flow_sparsifier_mwu.h"

#include "algorithms/semi_oblivious/expander_hierarchy/tree_sparsifier_solver.h"
#include "routing/storage/linear_routing_table.h"
#include "algorithms/semi_oblivious/expander_hierarchy/preprocessing/hierarchy_preprocessor.h"

Result<void> FlowSparsifier::run(LinearRoutingTable &table) {
    this->metrics.iteration_count = 10;

    _mwu_weights.resize(graph.getNumDirectedEdges());
    for(int e = 0; e < graph.getNumDirectedEdges(); e++) {
        _mwu_weights[e] = graph.getEdgeCapacity(e);
    }

    for (int t = 0; t<this->metrics.iteration_count ; t++) {
        double current_max_load = std::numeric_limits<double>::max();
        LinearRoutingTable current_table;
        current_table.init(graph);


        ElectrifiedExpanderHierarchySolver eeh_solver(this->graph, root);
        eeh_solver.init(_mwu_weights);
        auto xcut_res = eeh_solver.run(current_table);
        if (!xcut_res) {
            return getError(xcut_res);
        }

        std::vector<double> current_rload(graph.getNumDirectedEdges(), 0.0);
        for (int e = 0; e<graph.getNumDirectedEdges(); e++) {
            const auto& [head, tail] = graph.getEdgeEndpoints(e);
            if (head >= tail) continue;

            double& rload = current_rload[e];
            int rev_e = graph.getAntiEdge(e);

            for (int i{0};i<  current_table.src_ids[e].size(); i++) {
                rload += current_table.src_flows[e][i];
            }
            for (int i{0};i<  current_table.src_ids[rev_e].size(); i++) {
                rload += current_table.src_flows[rev_e][i];
            }
            rload /= graph.getEdgeCapacity(e);
            current_rload[rev_e] = rload;
        }

        _rel_loads = current_rload;

        tables.push_back(current_table);
        computeNewDistances();
    }


    return {};
}

Result<void> FlowSparsifier::scaleFlowDown(LinearRoutingTable &table) {
    if (tables.empty()) {
        return makeErrorMessage(ErrorCode::LogicError, "Computed zero tables.");
    }
    for (const auto& current_table : this->tables) {
        for (int e = 0; e<graph.getNumDirectedEdges(); e++) {
            for (int i{0};i<  current_table.src_ids[e].size(); i++) {
                table.addFlow(e, current_table.src_ids[e][i], current_table.src_flows[e][i]/this->metrics.iteration_count );
            }
        }
    }
    return {};
}
