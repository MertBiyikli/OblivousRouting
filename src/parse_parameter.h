//
// Created by Mert Biyikli on 23.08.25.
//

#ifndef OBLIVIOUSROUTING_PARSE_PARAMETER_H
#define OBLIVIOUSROUTING_PARSE_PARAMETER_H
#include <string>
#include <optional>
#include <algorithm>
#include "datastructures/GraphADJ.h"
#include "datastructures/GraphCSR.h"
#include "../experiments/performance/demands/DemandModel.h"
#include "../experiments/performance/demands/BimodalModel.h"
#include "../experiments/performance/demands/GaussianModel.h"
#include "../experiments/performance/demands/GravityModel.h"
#include "../experiments/performance/demands/UniformModel.h"
#include "lp_solver/MCCF_lp_solver.h"

enum class SolverType {
    ELECTRICAL_NAIVE,
    RAECKE_FRT,            // tree-based
    LP_APPLEGATE_COHEN,
    RAECKE_RANDOM_MST,
    ELECTRICAL_OPTIMIZED,
    ELECTRICAL_PARALLEL_BATCHES,
    ELECTRICAL_PARALLEL_ONTHEFLY,
    RAECKE_CKR,
    RAECKE_CKR_OPTIMIZED,
    OPTIMIZED_RAECKE_FRT,
};

enum class DemandModelType {
    GRAVITY,
    BIMODAL,
    GAUSSIAN,
    UNIFORM,
    NONE
};

struct Config {
    SolverType  solver;
    std::string filename;
    DemandModelType demand_model; // either "gravity" or "binomial" if a third argument is given
};

inline std::string to_lower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c){ return std::tolower(c); });
    return s;
}

inline std::optional<SolverType> parse_solver_token(std::string s) {
    s = to_lower(std::move(s));
    // electrical
    if (s == "electrical" || s == "electrical_naive" || s == "ef" || s == "e")
        return SolverType::ELECTRICAL_NAIVE;
    // tree / Räcke–FRT
    if (s == "tree" || s == "raecke" || s == "frt" || s == "r" || s == "t")
        return SolverType::RAECKE_FRT;
    // LP (Applegate–Cohen)
    if (s == "cohen" || s == "lp" || s == "applegate" || s == "ac" || s == "l")
        return SolverType::LP_APPLEGATE_COHEN;
    if (s == "mst" || s == "random_mst" || s == "raecke_mst" || s == "rmst")
        return SolverType::RAECKE_RANDOM_MST;
    if (s == "electrical_optimized" || s == "electricalopt" || s == "eo")
        return SolverType::ELECTRICAL_OPTIMIZED;
    if (s == "electrical_parallel_batches" || s == "electrical_batches" || s == "e_bacthes")
        return SolverType::ELECTRICAL_PARALLEL_BATCHES;
    if (s == "electrical_parallel_onthefly" || s == "electrical_onthefly" || s == "e_onthefly")
        return SolverType::ELECTRICAL_PARALLEL_ONTHEFLY;
    if (s == "ckr" || s == "ckr_partition" || s == "raecke_ckr")
        return SolverType::RAECKE_CKR;
    if (s == "ckr_optimized" || s == "ckropt" || s == "raecke_ckr_optimized")
        return SolverType::RAECKE_CKR_OPTIMIZED;
    if (s == "optimized_raecke_frt" || s == "raecke_frt_opt" || s == "oref")
        return SolverType::OPTIMIZED_RAECKE_FRT;


    if (s == "0") return SolverType::ELECTRICAL_NAIVE;
    if (s == "1") return SolverType::RAECKE_FRT;
    if (s == "2") return SolverType::LP_APPLEGATE_COHEN;
    if (s == "3") return SolverType::RAECKE_RANDOM_MST;
    if (s == "4") return SolverType::ELECTRICAL_OPTIMIZED;
    if (s == "5") return SolverType::ELECTRICAL_PARALLEL_BATCHES;
    if (s == "6") return SolverType::ELECTRICAL_PARALLEL_ONTHEFLY;
    if (s == "7") return SolverType::RAECKE_CKR;
    if (s == "8") return SolverType::RAECKE_CKR_OPTIMIZED;
    if (s == "9") return SolverType::OPTIMIZED_RAECKE_FRT;

    return std::nullopt;
}

inline std::optional<DemandModelType> parse_demand_model_token(std::string s) {
    s = to_lower(std::move(s));
    // gravity
    if (s == "gravity" || s == "gravity_model")
        return DemandModelType::GRAVITY;
    // binomial
    if (s == "bimodal" || s == "bimodal_model")
        return DemandModelType::BIMODAL;
    // gaussian
    if (s == "gaussian" || s == "gaussian_model")
        return DemandModelType::GAUSSIAN;
    // uniform
    if (s == "uniform" || s == "uniform_model")
        return DemandModelType::UNIFORM;

    return std::nullopt;
}

inline std::string usage(const char* prog) {
    std::ostringstream os;
    os << "Usage:\n"
       << "  " << prog << " <solver> <graph_file>\n\n"
       << "Solvers (case-insensitive):\n"
       << "  electrical | ef | e           -> Electrical Flow (naive)\n"
       << "  tree | raecke | frt | r | t   -> Tree-based (Raecke/FRT)\n"
       << "  cohen | lp | applegate | ac   -> Tree-based (Raecke/MST)\n"
       << "  mst | random_mst | raecke_mst | rmst -> LP (Raecke/Random MST)\n"
       << "  electrical_optimized | electricalopt | eo -> Electrical Flow (optimized)\n"
        << "  electrical_parallel_batches | electricalpar_batches | e_batches -> Electrical Flow (parallel)\n"
        << "  electrical_parallel_onthefly | electricalpar_onthefly | e_onthefly -> Electrical Flow (parallel)\n"
        << "  ckr | ckr_partition | raecke_ckr | -> Tree-based (Raecke/CKR)\n"
        << "  ckr_optimized | ckropt | raecke_ckr_optimized -> Tree-based (Raecke/CKR Optimized)\n\n"
       << "Numeric shortcuts: 0=electric, 1=tree, 2=cohen, 3=mst, 4=electrical_optimized_seq, 5=electrical_optimized_par_batches, 6=electrical_optimized_par_onthefly\n"
       << "[Optional] demand model (case-insensitive):\n"
       << "  gravity | gravity_model       -> Gravity Model\n"
       << "  binomial | binomial_model     -> Binomial Model\n"
       << "  gaussian | gaussian_model     -> Gaussian Model\n"
       << "  uniform | uniform_model       -> Uniform Model\n"
       << "If a demand model is provided, the oblivious ratio will be compared to the one of the demand model.\n";

    return os.str();
}

// Returns Config on success; prints an error to `err` string on failure.
inline std::optional<Config> parse_parameter(int argc, char** argv, std::string* err) {
    if (argc < 3) {
        if (err) *err = usage(argv[0]);
        return std::nullopt;
    }
    if (argc == 3) {
        auto solver_opt = parse_solver_token(argv[1]);
        if (!solver_opt) {
            if (err) *err = "Unknown solver: " + std::string(argv[1]) + "\n" + usage(argv[0]);
            return std::nullopt;
        }

        Config cfg{*solver_opt, std::string(argv[2]), DemandModelType::NONE}; // default demand model

        return cfg;
    }

    // if a third argument is given we use this to compare the solver to compare the oblivious ratio to
    // to one of the demand model, e.g. Gravity/Binomial model
    if ( argc == 4) {
        auto solver_opt = parse_solver_token(argv[1]);
        if (!solver_opt) {
            if (err) *err = "Unknown solver: " + std::string(argv[1]) + "\n" + usage(argv[0]);
            return std::nullopt;
        }

        auto demand_model_opt = parse_demand_model_token(argv[3]);
        if (!demand_model_opt) {
            if (err) *err = "Unknown demand model: " + std::string(argv[3]) + "\n" + usage(argv[0]);
            return std::nullopt;
        }

        Config cfg{*solver_opt, std::string(argv[2]), (demand_model_opt ? *demand_model_opt : DemandModelType::NONE)};
        return cfg;
    }else{
        if (err) *err = usage(argv[0]);
        return std::nullopt;
    }
}

inline std::pair<double, double> HandleDemandModel(int argc,
                              char** argv,
                              const std::optional<Config>& cfg,
                              GraphCSR& _g,
                              std::unique_ptr<ObliviousRoutingSolver>& _solver)
{
    if (argc < 4 || !cfg) { return {}; }

    std::cout << "Computing oblivious ratio for demand model: " << (argv[3] ? argv[3] : "<empty>") << std::endl;

    // if a demand model is provided, compute the oblivious ratio for that demand model
    if (cfg->demand_model != DemandModelType::NONE) {
        std::vector< std::pair<int, int> > demands;
        // print out the capacities for each node
        for (const auto& v : _g.getVertices()) {
            for (const auto& u : _g.getVertices()) {
                if (v != u) {
                    demands.push_back({v, u});
                }
            }
        }

        // initialize demand model
        std::unique_ptr<DemandModel> demand_model;
        switch  (cfg->demand_model) {
            case DemandModelType::GRAVITY:
                demand_model = std::make_unique<GravityModel>();
                break;
            case DemandModelType::BIMODAL:
                demand_model = std::make_unique<BimodalModel>();
                break;
            case DemandModelType::GAUSSIAN:
                demand_model = std::make_unique<GaussianModel>();
                break;
            case DemandModelType::UNIFORM:
                demand_model = std::make_unique<UniformModel>();
                break;
            default:
                std::cerr << "Unknown demand model type.\n";
        }
        // TODO: uncomment this
        /*
        std::cout << "Number of nodes: " << _g.getNumNodes() << std::endl;
        auto demand_map = demand_model->generate(_g, demands,  1.0);
        double total_demand = 0.0;
        for (const auto& [_, d] : demand_map) {
            total_demand += d;
        }

        // compute maximum congestion generated by demand
        std::vector congestion_per_edge(_g.getNumEdges()/2, 0.0); // undirected edges
        for (const auto& [edge, flow_map] : _solver->f_e_st) {
            double edge_cong = 0.0;
            for (const auto& [commodity, value] : flow_map) {
                edge_cong += std::abs(value * demand_map[commodity]);
            }

            edge_cong /= _g.getEdgeCapacity(edge.first, edge.second);
            int undirected_edge = INVALID_EDGE_ID;
            if (edge.first > edge.second) {
                // find the reverse edge
                undirected_edge = _g.getEdgeId(edge.second, edge.first);
            }else {
                undirected_edge = _g.getEdgeId(edge.first, edge.second);
            }
            if (undirected_edge != INVALID_EDGE_ID) {
                congestion_per_edge[undirected_edge] += edge_cong;
            }else {
                throw std::runtime_error("Anti edge not found in graph.");
            }
        }

        double max_cong = 0.0;
        for (const auto& cong : congestion_per_edge) {
            if (cong > max_cong) {
                max_cong = cong;
            }
        }


        std::cout << "Worst case congestion for " << (argv[3] ? argv[3] : "<empty>") << " model demand: " << max_cong << std::endl;



        // compute the Maximum Commodity flow for the given demand set
        CMMF_Solver mccf;
        mccf.setGraph(_g);
        mccf.init(_g);
        for (const auto& [d, value] : demand_map) {
            mccf.AddDemands({d.first, d.second}, value);
        }
        mccf.solve();
        double offline = mccf.getCongestion();

        return {offline, max_cong};
*/
    }
    return {};
}
#endif //OBLIVIOUSROUTING_PARSE_PARAMETER_H