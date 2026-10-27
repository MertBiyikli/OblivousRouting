//
// Created by Mert Biyikli on 23.06.26.
//

#ifndef OBLIVIOUSROUTING_TYPES_H
#define OBLIVIOUSROUTING_TYPES_H

#include <string>

enum class SolverType {
    ELECTRICAL_NAIVE,
    ELECTRICAL_SKETCHING,
    RAECKE_FRT_FLAT,
    RAECKE_CKR_FLAT,
    RAECKE_RANDOM_MST_FLAT,
    RAECKE_FRT_MENDELSCALING_FLAT,
    RAECKE_CKR_MENDELSCALING_FLAT,
    LP_APPLEGATE_COHEN,
    ELECTRICAL_PARALLEL_BATCHES,
    RAECKE_FRT_POINTER,
    RAECKE_CKR_POINTER,
    RAECKE_RANDOM_MST_POINTER,
    RAECKE_FRT_MENDELSCALING_POINTER,
    RAECKE_CKR_MENDELSCALING_POINTER,
    SEMI_ELECTRICAL,
    SEMI_TREE,
    EXPANDER_HIERARCHY,
    SEMI_EXPANDER_HIERARCHY,
    EXPANDER_MWU
};

enum class DemandModelType {
    GRAVITY,
    BIMODAL,
    GAUSSIAN,
    UNIFORM,
    NONE
};

enum class GraphFormat {
    CSR,
    ADJLIST
};


enum class OutputFormat {
    TEXT,
    JSON,
    COUT
};

enum class ResultStatus {
    OK,
    ERROR_INVALID_SOLVER,
    ERROR_INVALID_ROUTING_SCHEME,
    ERROR_MISSING_DEMAND_MODELS,
    ERROR_EXCEPTION
};


inline std::string getSolverTypeName(SolverType type) {
    switch (type) {
    case SolverType::ELECTRICAL_NAIVE:
        return "electrical_naive";

    case SolverType::ELECTRICAL_SKETCHING:
        return "electrical_sketching";

    case SolverType::RAECKE_FRT_FLAT:
        return "raecke_frt_flat";

    case SolverType::RAECKE_CKR_FLAT:
        return "raecke_ckr_flat";

    case SolverType::RAECKE_RANDOM_MST_FLAT:
        return "raecke_random_mst_flat";

    case SolverType::RAECKE_FRT_MENDELSCALING_FLAT:
        return "raecke_frt_mendel_flat";

    case SolverType::RAECKE_CKR_MENDELSCALING_FLAT:
        return "raecke_ckr_mendel_flat";

    case SolverType::LP_APPLEGATE_COHEN:
        return "lp_applegate_cohen";

    case SolverType::ELECTRICAL_PARALLEL_BATCHES:
        return "electrical_parallel";

    case SolverType::RAECKE_FRT_POINTER:
        return "raecke_frt_pointer";

    case SolverType::RAECKE_CKR_POINTER:
        return "raecke_ckr_pointer";

    case SolverType::RAECKE_RANDOM_MST_POINTER:
        return "raecke_random_mst_pointer";

    case SolverType::RAECKE_FRT_MENDELSCALING_POINTER:
        return "raecke_frt_mendel_pointer";

    case SolverType::RAECKE_CKR_MENDELSCALING_POINTER:
        return "raecke_ckr_mendel_pointer";

    case SolverType::SEMI_ELECTRICAL:
        return "semi_electrical";

    case SolverType::SEMI_TREE:
        return "semi_tree";

    case SolverType::EXPANDER_HIERARCHY:
        return "expander_hierarchy";

    case SolverType::SEMI_EXPANDER_HIERARCHY:
        return "semi_expander_hierarchy";

    case SolverType::EXPANDER_MWU:
        return "expander_mwu";
    }

    return "unknown";
}

#endif //OBLIVIOUSROUTING_TYPES_H