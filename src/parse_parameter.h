//
// Created by Mert Biyikli on 23.08.25.
//

#ifndef OBLIVIOUSROUTING_PARSE_PARAMETER_H
#define OBLIVIOUSROUTING_PARSE_PARAMETER_H
#include <string>
#include <optional>
#include <algorithm>
#include "../experiments/performance/demands/DemandModel.h"
#include "../experiments/performance/demands/BimodalModel.h"
#include "../experiments/performance/demands/GaussianModel.h"
#include "../experiments/performance/demands/GravityModel.h"
#include "../experiments/performance/demands/UniformModel.h"
#include "lp_solver/LPSolver.h"
#include "lp_solver/MCCF_lp_solver.h"

#include "electrical/electrical.h"
#include "electrical/parallel/electrical_flow_par_batches.h"

#include "tree_based/tree_mwu.h"
#include "tree_based/frt.h"
#include "tree_based/fast_ckr.h"
#include "tree_based/random_mst.h"

std::vector<std::string> solverNames ={
    "electrical",
    "raecke_frt",
    "raecke_ckr",
       "raecke_random_mst",
    "raecke_frt_mendel",
    "raecke_ckr_mendel",
    "lp_applegate_cohen",
    "electrical_parallel_batches"
};

enum class SolverType {
    ELECTRICAL,
    RAECKE_FRT,
    RAECKE_CKR,
    RAECKE_RANDOM_MST,
    RAECKE_FRT_MENDELSCALING,
    RAECKE_CKR_MENDELSCALING,
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


enum class GraphFormat {
    CSR,
    ADJLIST,
    ADJMATRIX
};

struct Config {
    std::vector<SolverType>  solvers;
    std::string filename;
    DemandModelType demand_model;
    GraphFormat graph_format;
};

std::unique_ptr<ObliviousRoutingSolver>
makeSolver(SolverType type, IGraph& g) {
    switch (type) {
        case SolverType::ELECTRICAL:
            return std::make_unique<ElectricalMWU>(g, 0);

        case SolverType::RAECKE_FRT:
            return std::make_unique<TreeMWU>(g, 0, std::make_unique<FRT>(g));

        case SolverType::RAECKE_CKR:
            return std::make_unique<TreeMWU>(g, 0, std::make_unique<FastCKR>(g));

        case SolverType::RAECKE_RANDOM_MST:
            return std::make_unique<TreeMWU>(g,0,  std::make_unique<TreeMST>(g));

        case SolverType::RAECKE_FRT_MENDELSCALING:
            return std::make_unique<TreeMWU>(g, 0, std::make_unique<FRT>(g, true));

        case SolverType::RAECKE_CKR_MENDELSCALING:
            return std::make_unique<TreeMWU>(g, 0, std::make_unique<FastCKR>(g, true));

        case SolverType::LP_APPLEGATE_COHEN:
            return std::make_unique<LPSolver>(g);

        case SolverType::ELECTRICAL_PARALLEL_BATCHES:
            return std::make_unique<ParallelElectricalMWU>(g, 0);

        default:
            throw std::runtime_error("Unknown solver type.");
    }
}


std::unique_ptr<IGraph> makegraph(GraphFormat type) {
    switch (type) {
        case GraphFormat::CSR:
            return std::make_unique<GraphCSR>();

        case GraphFormat::ADJLIST:
            return std::make_unique<GraphADJList>();

        case GraphFormat::ADJMATRIX:
            return std::make_unique<GraphADJMatrix>();

        default:
            throw std::runtime_error("Unknown graph format.");
    }
}


inline std::string to_lower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c){ return std::tolower(c); });
    return s;
}

inline std::optional<SolverType> parse_solver_token(std::string s) {
    s = to_lower(std::move(s));

    if (s == "electrical" || s == "elec" || s == "ef" || s == "e")
        return SolverType::ELECTRICAL;
    if (s == "raecke_frt" || s == "frt" || s == "f")
        return SolverType::RAECKE_FRT;
    if (s == "raecke_ckr" || s == "ckr" || s == "c")
        return SolverType::RAECKE_CKR;
    if (s == "raecke_mst" || s == "random_mst" || s == "rmst" || s == "mst")
        return SolverType::RAECKE_RANDOM_MST;
    if (s == "raecke_frt_mendel" || s == "frt_mendel")
        return SolverType::RAECKE_FRT_MENDELSCALING;
    if (s == "raecke_ckr_mendel" || s == "ckr_mendel")
        return SolverType::RAECKE_CKR_MENDELSCALING;
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
    if (s == "6") return SolverType::RAECKE_FRT_MENDELSCALING;
    if (s == "7") return SolverType::RAECKE_CKR_MENDELSCALING;

    return std::nullopt;
}

inline std::optional<DemandModelType> parse_demand_model_token(std::string s) {
    s = to_lower(std::move(s));
    if (s == "gravity" || s == "gravity_model")
        return DemandModelType::GRAVITY;
    if (s == "bimodal" || s == "bimodal_model")
        return DemandModelType::BIMODAL;
    if (s == "gaussian" || s == "gaussian_model")
        return DemandModelType::GAUSSIAN;
    if (s == "uniform" || s == "uniform_model")
        return DemandModelType::UNIFORM;

    return std::nullopt;
}


inline std::optional<GraphFormat> parse_graph_format_token(std::string s) {
    s = to_lower(std::move(s));
    if (s == "csr") return GraphFormat::CSR;
    if (s == "adjlist" || s == "list") return GraphFormat::ADJLIST;
    if (s == "adjmatrix" || s == "matrix") return GraphFormat::ADJMATRIX;

    return std::nullopt;
}

inline std::string usage(const char* prog) {
    std::ostringstream os;
    os << "Usage:\n"
       << "  " << "./oblivious_routing " << " <solver> <graph_file> OPTIOnAL: <demand_model> <graph_format> \n\n"
       << "Solvers (case-insensitive):\n"
       << "  electrical | ef | e           -> Electrical Flow (naive)\n"
       << "  raecke_frt | frt | f   -> Tree-based (Raecke/FRT)\n"
        << " raecke_ckr | ckr | c -> Tree-based (Raecke/CKR)\n"
        << "  raecke_frt_mendel | frt_mendel   -> Tree-based (Raecke/FRT) using MendelScaling\n"
        << " raecke_ckr_mendel | ckr_mendel -> Tree-based (Raecke/CKR) using MendelScaling\n"
       << "  mst | random_mst | raecke_mst | rmst -> LP (Raecke/Random MST)\n"
       << "  cohen | lp | applegate | ac   -> Tree-LP (LP/Applegate and Cohen)\n"
        << "  electrical_parallel | elec_par | e_par -> Electrical Flow (parallel)\n"
       << "Numeric shortcuts: 0=electric, 1=frt, 2=ckr, 3=mst, 4=LP, 5=electrical_parallel, 6=frt_mendel, 7=ckr_mendel\n"
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
            DemandModelType::NONE,
            GraphFormat::CSR
        };
        return cfg;
    }

    if (argc == 4) {
        auto demand_model_opt = parse_demand_model_token(argv[3]);
        if (!demand_model_opt) {
            auto graph_format_opt = parse_graph_format_token(argv[3]);
            if (graph_format_opt) {
                Config cfg{
                    *solvers_opt,
                    std::string(argv[2]),
                    DemandModelType::NONE,
                    *graph_format_opt
                };
                return cfg;
            }else {
                if (err) *err = "Unknown demand model/Graph format model: " + std::string(argv[3]) + "\n" + usage(argv[0]);
                return std::nullopt;
            }
        }

        Config cfg{
            *solvers_opt,
            std::string(argv[2]),
            *demand_model_opt,
            GraphFormat::CSR
        };
        return cfg;
    }

    if (argc == 5) {
        auto demand_model_opt = parse_demand_model_token(argv[3]);
        if (!demand_model_opt) {
            if (err) *err = "Unknown demand model: " + std::string(argv[3]) + "\n" + usage(argv[0]);
            return std::nullopt;
        }

        auto graph_format_opt = parse_graph_format_token(argv[4]);
        if (!graph_format_opt) {
            if (err) *err = "Unknown graph format: " + std::string(argv[4]) + "\n" + usage(argv[0]);
            return std::nullopt;
        }

        Config cfg{
            *solvers_opt,
            std::string(argv[2]),
            *demand_model_opt,
            *graph_format_opt
        };
        return cfg;
    }

    if (err) *err = usage(argv[0]);
    return std::nullopt;
}


inline std::pair<double, double> HandleDemandModel(int argc,
                              char** argv,
                              const std::optional<Config>& cfg,
                              IGraph& _g,
                              const std::unique_ptr<RoutingScheme>& routing_scheme)
{
    if (argc < 4 || !cfg) { return {}; }

    // find demand model type from config and print it
    std::string demand_model_str = (argv[3] ? argv[3] : "<empty>");
    auto demand_model = parse_demand_model_token(demand_model_str);
    if (demand_model) {
        std::cout << "Demand model: " << (argv[3] ? argv[3] : "<empty>") << std::endl;
    }
    // if a demand model is provided, compute the oblivious ratio for that demand model
    if (cfg->demand_model != DemandModelType::NONE) {
        std::vector< std::pair<int, int> > demands;

        for (const auto& v : _g.getVertices()) {
            for (const auto& u : _g.getVertices()) {
                if (v != u) {
                    demands.push_back({v, u});
                }
            }
        }

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

        std::vector<double> congestion_per_edge(_g.getNumEdges(), 0.0); // undirected edges

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


void printStatsForDemandModel(char** argv,std::pair<double, double> result) {
    if (result.first > 0.0 && result.second > 0.0) {
        std::cout << "Ratio off the optimal offline solution "
                  << (argv[3] ? argv[3] : "<empty>")
                  << " model demand: "
                  << ( result.second  / result.first ) * 100.0 << "% "
                  << " ("
                  << result.first << " / " << result.second
                  << ")\n";
    }else {
        std::cout << "Invalid congestion values for demand model evaluation.\n";
    }
}

#endif //OBLIVIOUSROUTING_PARSE_PARAMETER_H