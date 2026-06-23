//
// Created by Mert Biyikli on 25.03.26.
//

#ifndef OBLIVIOUSROUTING_SOLVER_IO_H
#define OBLIVIOUSROUTING_SOLVER_IO_H

#include "core/solver.h"
#include "algorithms/oblivious/oblivious_solver.h"
#include "algorithms/lp/lp_ac.h"
#include "algorithms/oblivious/mwu/electrical_mwu.h"
#include "algorithms/oblivious/mwu/tree_mwu.h"
#include "algorithms/oblivious/mwu/oracle/tree/mst/mst_oracle.h"
#include "algorithms/oblivious/mwu/oracle/tree/frt/frt.h"
#include "algorithms/oblivious/mwu/oracle/tree/fast_ckr/fast_ckr.h"

#include <string>
#include <optional>
#include <vector>
#include <memory>
#include <functional>
#include <map>

#include "algorithms/semi_oblivious/or_tools_optimizer.h"
#include "algorithms/semi_oblivious/semi_oblivious_solver.h"
#include "core/types.h"

// Map-based token parsers for reduced code duplication
static const std::map<std::string, SolverType> SOLVER_MAP{
    {"electrical_naive", SolverType::ELECTRICAL_NAIVE}, {"elec_naive", SolverType::ELECTRICAL_NAIVE},
    {"ef_naive", SolverType::ELECTRICAL_NAIVE},  {"0", SolverType::ELECTRICAL_NAIVE},
     {"electrical", SolverType::ELECTRICAL_SKETCHING}, {"electrical_sketching", SolverType::ELECTRICAL_SKETCHING}, {"elec_sketching", SolverType::ELECTRICAL_SKETCHING},
    {"ef_sketching", SolverType::ELECTRICAL_SKETCHING}, {"e_sketching", SolverType::ELECTRICAL_SKETCHING}, {"1", SolverType::ELECTRICAL_SKETCHING},
    {"raecke_frt", SolverType::RAECKE_FRT_FLAT}, {"frt", SolverType::RAECKE_FRT_FLAT},
    {"f", SolverType::RAECKE_FRT_FLAT}, {"2", SolverType::RAECKE_FRT_FLAT},
    {"raecke_ckr", SolverType::RAECKE_CKR_FLAT}, {"ckr", SolverType::RAECKE_CKR_FLAT},
    {"c", SolverType::RAECKE_CKR_FLAT}, {"3", SolverType::RAECKE_CKR_FLAT},
    {"raecke_mst", SolverType::RAECKE_RANDOM_MST_FLAT}, {"random_mst", SolverType::RAECKE_RANDOM_MST_FLAT},
    {"rmst", SolverType::RAECKE_RANDOM_MST_FLAT}, {"mst", SolverType::RAECKE_RANDOM_MST_FLAT}, {"4", SolverType::RAECKE_RANDOM_MST_FLAT},
    {"cohen", SolverType::LP_APPLEGATE_COHEN}, {"lp", SolverType::LP_APPLEGATE_COHEN},
    {"applegate", SolverType::LP_APPLEGATE_COHEN}, {"ac", SolverType::LP_APPLEGATE_COHEN}, {"l", SolverType::LP_APPLEGATE_COHEN}, {"5", SolverType::LP_APPLEGATE_COHEN},
    {"electrical_parallel", SolverType::ELECTRICAL_PARALLEL_BATCHES}, {"elec_par", SolverType::ELECTRICAL_PARALLEL_BATCHES},
    {"e_par", SolverType::ELECTRICAL_PARALLEL_BATCHES}, {"6", SolverType::ELECTRICAL_PARALLEL_BATCHES},
    {"raecke_frt_mendel", SolverType::RAECKE_FRT_MENDELSCALING_FLAT}, {"frt_mendel", SolverType::RAECKE_FRT_MENDELSCALING_FLAT}, {"7", SolverType::RAECKE_FRT_MENDELSCALING_FLAT},
    {"raecke_ckr_mendel", SolverType::RAECKE_CKR_MENDELSCALING_FLAT}, {"ckr_mendel", SolverType::RAECKE_CKR_MENDELSCALING_FLAT}, {"8", SolverType::RAECKE_CKR_MENDELSCALING_FLAT},
    {"raecke_frt_pointer", SolverType::RAECKE_FRT_POINTER}, {"frt_pointer", SolverType::RAECKE_FRT_POINTER}, {"9", SolverType::RAECKE_FRT_POINTER},
    {"raecke_ckr_pointer", SolverType::RAECKE_CKR_POINTER}, {"ckr_pointer", SolverType::RAECKE_CKR_POINTER}, {"10", SolverType::RAECKE_CKR_POINTER},
    {"raecke_mst_pointer", SolverType::RAECKE_RANDOM_MST_POINTER}, {"random_mst_pointer", SolverType::RAECKE_RANDOM_MST_POINTER},
    {"rmst_pointer", SolverType::RAECKE_RANDOM_MST_POINTER}, {"mst_pointer", SolverType::RAECKE_RANDOM_MST_POINTER}, {"11", SolverType::RAECKE_RANDOM_MST_POINTER},
    {"raecke_frt_mendel_pointer", SolverType::RAECKE_FRT_MENDELSCALING_POINTER}, {"frt_mendel_pointer", SolverType::RAECKE_FRT_MENDELSCALING_POINTER}, {"12", SolverType::RAECKE_FRT_MENDELSCALING_POINTER},
        {"semi_elec", SolverType::SEMI_ELECTRICAL}, {"semi_electrical", SolverType::SEMI_ELECTRICAL}, {"13", SolverType::SEMI_ELECTRICAL},
        {"semi_tree", SolverType::SEMI_TREE}, {"14", SolverType::SEMI_TREE}
};


inline std::optional<std::unique_ptr<ISolver>>
makeSolver(SolverType type, IGraph& g) {
    // Factory with cycle removal strategy support for TreeMWU-based solvers
    switch (type) {
        case SolverType::ELECTRICAL_NAIVE:
            return std::make_unique<ElectricalMWU>(g, 0, false);

            case SolverType::ELECTRICAL_SKETCHING:
            return std::make_unique<ElectricalMWU>(g, 0, true);

        case SolverType::RAECKE_FRT_FLAT:
            return std::make_unique<TreeMWU<FlatHST>>(g, 0, std::make_unique<FRT<FlatHST>>(g, false));

        case SolverType::RAECKE_CKR_FLAT:
            return std::make_unique<TreeMWU<FlatHST>>(g, 0, std::make_unique<FastCKR<FlatHST>>(g, false));

        case SolverType::RAECKE_RANDOM_MST_FLAT:
            return std::make_unique<TreeMWU<FlatHST>>(g, 0, std::make_unique<TreeMST<FlatHST>>(g));

        case SolverType::RAECKE_FRT_MENDELSCALING_FLAT:
            return std::make_unique<TreeMWU<FlatHST>>(g, 0, std::make_unique<FRT<FlatHST>>(g, true));

        case SolverType::RAECKE_CKR_MENDELSCALING_FLAT:
            return std::make_unique<TreeMWU<FlatHST>>(g, 0, std::make_unique<FastCKR<FlatHST>>(g, true));

        case SolverType::LP_APPLEGATE_COHEN:
            return std::make_unique<LPSolver>(g);

        case SolverType::RAECKE_FRT_POINTER:
            return std::make_unique<TreeMWU<std::shared_ptr<HSTNode>>>(g, 0, std::make_unique<FRT<std::shared_ptr<HSTNode>>>(g, false));

        case SolverType::RAECKE_CKR_POINTER:
            return std::make_unique<TreeMWU<std::shared_ptr<HSTNode>>>(g, 0, std::make_unique<FastCKR<std::shared_ptr<HSTNode>>>(g, false));

        case SolverType::RAECKE_RANDOM_MST_POINTER:
            return std::make_unique<TreeMWU<std::shared_ptr<HSTNode>>>(g, 0, std::make_unique<TreeMST<std::shared_ptr<HSTNode>>>(g));

        case SolverType::RAECKE_FRT_MENDELSCALING_POINTER:
            return std::make_unique<TreeMWU<std::shared_ptr<HSTNode>>>(g, 0, std::make_unique<FRT<std::shared_ptr<HSTNode>>>(g, true));

        case SolverType::RAECKE_CKR_MENDELSCALING_POINTER:
            return std::make_unique<TreeMWU<std::shared_ptr<HSTNode>>>(g, 0, std::make_unique<FastCKR<std::shared_ptr<HSTNode>>>(g, true));

        default:
            return std::nullopt;
    }
}




// Helper to get solver name from type
inline std::string getSolverName(SolverType type) {
    static const std::map<SolverType, std::string> names{
            {SolverType::ELECTRICAL_NAIVE, "Electrical Flow (naive)"},
            {SolverType::ELECTRICAL_SKETCHING, "Electrical Flow (sketching)"},
            {SolverType::RAECKE_FRT_FLAT, "Raecke FRT (Flat HST)"},
            {SolverType::RAECKE_CKR_FLAT, "Raecke CKR (Flat HST)"},
            {SolverType::RAECKE_RANDOM_MST_FLAT, "Random MST (Flat HST)"},
            {SolverType::RAECKE_FRT_MENDELSCALING_FLAT, "Raecke FRT + MendelScaling (Flat HST)"},
            {SolverType::RAECKE_CKR_MENDELSCALING_FLAT, "Raecke CKR + MendelScaling (Flat HST)"},
            {SolverType::LP_APPLEGATE_COHEN, "LP Applegate-Cohen"},
            {SolverType::ELECTRICAL_PARALLEL_BATCHES, "Electrical Flow (Parallel)"},
            {SolverType::RAECKE_FRT_POINTER, "Raecke FRT (Pointer HST)"},
            {SolverType::RAECKE_CKR_POINTER, "Raecke CKR (Pointer HST)"},
            {SolverType::RAECKE_RANDOM_MST_POINTER, "Random MST (Pointer HST)"},
            {SolverType::RAECKE_FRT_MENDELSCALING_POINTER, "Raecke FRT + MendelScaling (Pointer HST)"},
            {SolverType::RAECKE_CKR_MENDELSCALING_POINTER, "Raecke CKR + MendelScaling (Pointer HST)"},
            {SolverType::SEMI_ELECTRICAL, "Semi-Oblivious Routing (Electrical base)"},
            {SolverType::SEMI_TREE, "Semi-Oblivious Routing (Tree base)"}
    };
    auto it = names.find(type);
    return (it != names.end()) ? it->second : "Unknown Solver";
}



#endif //OBLIVIOUSROUTING_SOLVER_IO_H