//
// Created by Mert on 07.10.26.
//

#ifndef E_ROUTING_FAILURE_RECOVER_ANALYSIS_H
#define E_ROUTING_FAILURE_RECOVER_ANALYSIS_H

#include <cstddef>
#include <vector>

#include "core/errors.h"
#include "core/types.h"
#include "data_structures/graph/graph.h"
#include "routing/routing_table.h"
#include "utils/demands.h"
#include "algorithms/semi_oblivious/semi_oblivious_solver.h"
#include "algorithms/semi_oblivious/semi_routing_engine.h"
#include "algorithms/semi_oblivious/postprocessing/or_tools_optimizer.h"
#include "visualization/failure_recovery_types.h"


class FailureRecoveryAnalyzer {
public:

    /*
     * Performs complete N-1 physical-link recovery analysis.
     *
     * Currently intended for normal oblivious solvers:
     *
     *   ELECTRICAL_SKETCHING
     *   RAECKE_CKR_FLAT
     *   RAECKE_FRT_FLAT
     *   EXPANDER_HIERARCHY
     *   ...
     *
     * Semi-oblivious solvers will be added after this path is
     * validated because their recovery lifecycle also includes
     * candidate-path preprocessing + demand-specific optimization.
     */
    static Result<FailureRecoveryAnalysis> analyze(optimized::Graph<EdgeData>& graph, SolverType solver_type, const demands& demand_map, const DemandModelType& demand_model, double baseline_congestion);

private:
    static Result<std::unique_ptr<RoutingScheme>> recomputeSemiOblivious(
        optimized::Graph<EdgeData>& graph,
        SolverType solver_type,
        const demands& demand_map,
        DemandModelType demand_model
    );

    static bool isSemiObliviousSolver(SolverType solver_type);
};

#endif //E_ROUTING_FAILURE_RECOVER_ANALYSIS_H
