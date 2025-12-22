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
#include "lp_solver/LPSolver.h"
#include "lp_solver/MCCF_lp_solver.h"
#include "tree_based/ckr/raecke_mwu_ckr.h"
#include "tree_based/frt/raecke_mwu_frt.h"
#include "tree_based/random_mst/raecke_mwu_random.h"

std::vector<std::string> solverNames ={
    "electrical",
    "raecke_frt",
    "raecke_ckr",
    "raecke_random_mst",
    "lp_applegate_cohen",
    "electrical_parallel_batches"
};

enum class SolverType {
    ELECTRICAL,
    RAECKE_FRT,            // tree-based
    RAECKE_CKR,
    RAECKE_RANDOM_MST,
    LP_APPLEGATE_COHEN,
    ELECTRICAL_PARALLEL_BATCHES
};

enum class DemandModelType {
    GRAVITY,
    BIMODAL,
    GAUSSIAN,
    UNIFORM,
    NONE
};

struct Config {
    std::vector<SolverType>  solvers;
    std::string filename;
    DemandModelType demand_model; // either "gravity" or "binomial" if a third argument is given
};

std::unique_ptr<ObliviousRoutingSolver>
makeSolver(SolverType type, IGraph& g_csr) {
    switch (type) {
        case SolverType::ELECTRICAL:
            return std::make_unique<ElectricalFlowOptimized>(g_csr, 0);

        case SolverType::RAECKE_FRT:
            return std::make_unique<RaeckeMWU_FRT>(g_csr);

        case SolverType::RAECKE_CKR:
            return std::make_unique<RaeckeMWU_CKR>(g_csr);

        case SolverType::RAECKE_RANDOM_MST:
            return std::make_unique<RaeckeMWU_Random>(g_csr);

        case SolverType::LP_APPLEGATE_COHEN:
            return std::make_unique<LPSolver>(g_csr);

        case SolverType::ELECTRICAL_PARALLEL_BATCHES:
            // return std::make_unique<ElectricalFlowParallelBatches>(g_csr);
            throw std::runtime_error("Parallel batches solver not implemented.");

        default:
            throw std::runtime_error("Unknown solver type.");
    }
}


inline std::string to_lower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c){ return std::tolower(c); });
    return s;
}

inline std::optional<SolverType> parse_solver_token(std::string s) {
    s = to_lower(std::move(s));
    // electrical
    if (s == "electrical" || s == "elec" || s == "ef" || s == "e")
        return SolverType::ELECTRICAL;
    // tree / Räcke–FRT
    if (s == "raecke_frt" || s == "frt" || s == "f")
        return SolverType::RAECKE_FRT;
    if (s == "raecke_ckr" || s == "ckr" || s == "c")
        return SolverType::RAECKE_CKR;
    if (s == "raecke_mst" || s == "random_mst" || s == "rmst")
        return SolverType::RAECKE_RANDOM_MST;
    // LP (Applegate–Cohen)
    if (s == "cohen" || s == "lp" || s == "applegate" || s == "ac" || s == "l")
        return SolverType::LP_APPLEGATE_COHEN;
    if (s == "electrical_parallel" || s == "elec_par" || s == "e_par")
        return SolverType::ELECTRICAL_PARALLEL_BATCHES;



    if (s == "0") return SolverType::ELECTRICAL;
    if (s == "1") return SolverType::RAECKE_FRT;
    if (s == "2") return SolverType::RAECKE_CKR;
    if (s == "3") return SolverType::RAECKE_RANDOM_MST;
    if (s == "4") return SolverType::LP_APPLEGATE_COHEN;
    if (s == "5") return SolverType::ELECTRICAL_PARALLEL_BATCHES;

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
       << "  raecke_frt | frt | f   -> Tree-based (Raecke/FRT)\n"
        << " raecke_ckr | ckr | c -> Tree-based (Raecke/CKR)\n"
       << "  cohen | lp | applegate | ac   -> Tree-based (LP/Applegate and Cohen)\n"
       << "  mst | random_mst | raecke_mst | rmst -> LP (Raecke/Random MST)\n"
        << "  electrical_parallel | elec_par | e_par -> Electrical Flow (parallel)\n"
       << "Numeric shortcuts: 0=electric, 1=frt, 2=ckr, 3=mst, 4=LP, 5=electrical_parallel\n"
       << "[Optional] demand model (case-insensitive):\n"
       << "  gravity | gravity_model       -> Gravity Model\n"
       << "  binomial | binomial_model     -> Binomial Model\n"
       << "  gaussian | gaussian_model     -> Gaussian Model\n"
       << "  uniform | uniform_model       -> Uniform Model\n"
       << "If a demand model is provided, the oblivious ratio will be compared to the one of the demand model.\n";

    return os.str();
}


inline std::optional<std::vector<SolverType>>
parse_solver_list(std::string s) {
    std::vector<SolverType> result;

    size_t start = 0;
    while (true) {
        size_t pos = s.find(',', start);
        std::string token = (pos == std::string::npos)
                                ? s.substr(start)
                                : s.substr(start, pos - start);

        auto solver = parse_solver_token(token);
        if (!solver) return std::nullopt;

        result.push_back(*solver);

        if (pos == std::string::npos) break;
        start = pos + 1;
    }

    if (result.empty()) return std::nullopt;
    return result;
}


// Returns Config on success; prints an error to `err` string on failure.
inline std::optional<Config> parse_parameter(int argc, char** argv, std::string* err) {
    if (argc < 3) {
        if (err) *err = usage(argv[0]);
        return std::nullopt;
    }

    auto solvers_opt = parse_solver_list(argv[1]);
    if (!solvers_opt) {
        if (err) *err = "Unknown solver list: " + std::string(argv[1]) + "\n" + usage(argv[0]);
        return std::nullopt;
    }

    if (argc == 3) {
        Config cfg{
            *solvers_opt,
            std::string(argv[2]),
            DemandModelType::NONE
        };
        return cfg;
    }

    if (argc == 4) {
        auto demand_model_opt = parse_demand_model_token(argv[3]);
        if (!demand_model_opt) {
            if (err) *err = "Unknown demand model: " + std::string(argv[3]) + "\n" + usage(argv[0]);
            return std::nullopt;
        }

        Config cfg{
            *solvers_opt,
            std::string(argv[2]),
            *demand_model_opt
        };
        return cfg;
    }

    if (err) *err = usage(argv[0]);
    return std::nullopt;
}


inline std::pair<double, double> HandleDemandModel(int argc,
                              char** argv,
                              const std::optional<Config>& cfg,
                              GraphCSR& _g,
                              const std::unique_ptr<RoutingScheme>& routing_scheme)
{
    if (argc < 4 || !cfg) { return {}; }

    std::cout << "Demand model: " << (argv[3] ? argv[3] : "<empty>") << std::endl;

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

        auto demand_map = demand_model->generate(_g, demands,  1.0);
        double total_demand = 0.0;
        for (int i = 0; i < demand_map.size(); ++i) {
            total_demand += demand_map.getDemandValue(i);
        }

        // compute maximum congestion generated by demand
        std::vector congestion_per_edge(_g.getNumEdges(), 0.0); // undirected edges

        routing_scheme->routeDemands(congestion_per_edge, demand_map);
        double max_cong = routing_scheme->getMaxCongestion(congestion_per_edge);



        // compute the Maximum Commodity flow for the given demand set
        CMMF_Solver mccf(_g);
        mccf.init(_g);
        mccf.AddDemandMap(demand_map);
        auto offline_scheme = mccf.solve();

        double offline_congestion = mccf.getCongestion();
        return {offline_congestion, max_cong};

    }
    return {};
}
#endif //OBLIVIOUSROUTING_PARSE_PARAMETER_H