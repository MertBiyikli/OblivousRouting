//
// Created by Mert Biyikli on 25.03.26.
//

#ifndef OBLIVIOUSROUTING_SOLVER_IO_H
#define OBLIVIOUSROUTING_SOLVER_IO_H

#include "../include/algorithms/lp/lp_ac.h"
#include "../include/algorithms/mwu/electrical_mwu.h"
#include "../include/algorithms/mwu/tree_mwu.h"
#include "../include/algorithms/mwu/oracle/tree/mst/mst_oracle.h"
#include "../include/algorithms/mwu/oracle/tree/frt/frt.h"
#include "../include/algorithms/mwu/oracle/tree/fast_ckr/fast_ckr.h"

#include <string>
#include <optional>
#include <vector>
#include <memory>
#include <functional>
#include <map>

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
    RAECKE_CKR_MENDELSCALING_POINTER
};

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
};


inline std::optional<std::unique_ptr<ObliviousRoutingSolver>>
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

#endif //OBLIVIOUSROUTING_SOLVER_IO_H