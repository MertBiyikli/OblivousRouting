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
#include "tree_based/utils/quotient_graph.h"
#include "tree_based/utils/ultrametric_tree.h"

std::vector<std::string> solverNames ={
    "electrical",
    "raecke_frt",
    "raecke_ckr",
    "raecke_random_mst",
    "raecke_frt_mendel",
    "raecke_ckr_mendel",
    "lp_applegate_cohen",
    "electrical_parallel_batches",
    "racke_frt_pointer",
    "racke_ckr_pointer",
    "racke_random_mst_pointer",
    "racke_frt_mendel_pointer",
    "racke_ckr_mendel_pointer"
};

enum class SolverType {
    ELECTRICAL,
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
    std::vector<SolverType>      solvers;
    std::string                  filename;
    std::vector<DemandModelType> demand_models;  // empty = no demand evaluation
    GraphFormat                  graph_format;
};

std::unique_ptr<ObliviousRoutingSolver>
makeSolver(SolverType type, IGraph& g) {
    switch (type) {
        case SolverType::ELECTRICAL:
            return std::make_unique<ElectricalMWU>(g, 0);

        case SolverType::RAECKE_FRT_FLAT:
            return std::make_unique<TreeMWU<FlatHST>>(g, 0, std::make_unique<FRT<FlatHST>>(g,false));

        case SolverType::RAECKE_CKR_FLAT:
            return std::make_unique<TreeMWU<FlatHST>>(g, 0, std::make_unique<FastCKR<FlatHST>>(g,false));

        case SolverType::RAECKE_RANDOM_MST_FLAT:
            return std::make_unique<TreeMWU<FlatHST>>(g,0,  std::make_unique<TreeMST<FlatHST>>(g));

        case SolverType::RAECKE_FRT_MENDELSCALING_FLAT:
            return std::make_unique<TreeMWU<FlatHST>>(g, 0, std::make_unique<FRT<FlatHST>>(g, true));

        case SolverType::RAECKE_CKR_MENDELSCALING_FLAT:
            return std::make_unique<TreeMWU<FlatHST>>(g, 0, std::make_unique<FastCKR<FlatHST>>(g, true));

        case SolverType::LP_APPLEGATE_COHEN:
            return std::make_unique<LPSolver>(g);

        case SolverType::ELECTRICAL_PARALLEL_BATCHES:
            return std::make_unique<ParallelElectricalMWU>(g, 0);

        case SolverType::RAECKE_FRT_POINTER:
            return std::make_unique<TreeMWU<std::shared_ptr<HSTNode>>>(g, 0, std::make_unique<FRT<std::shared_ptr<HSTNode>>>(g,false));

        case SolverType::RAECKE_CKR_POINTER:
            return std::make_unique<TreeMWU<std::shared_ptr<HSTNode>>>(g, 0, std::make_unique<FastCKR<std::shared_ptr<HSTNode>>>(g,false));

        case SolverType::RAECKE_RANDOM_MST_POINTER:
            return std::make_unique<TreeMWU<std::shared_ptr<HSTNode>>>(g,0,  std::make_unique<TreeMST<std::shared_ptr<HSTNode>>>(g));

        case SolverType::RAECKE_FRT_MENDELSCALING_POINTER:
            return std::make_unique<TreeMWU<std::shared_ptr<HSTNode>>>(g, 0, std::make_unique<FRT<std::shared_ptr<HSTNode>>>(g, true));

        case SolverType::RAECKE_CKR_MENDELSCALING_POINTER:
            return std::make_unique<TreeMWU<std::shared_ptr<HSTNode>>>(g, 0, std::make_unique<FastCKR<std::shared_ptr<HSTNode>>>(g, true));


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
        return SolverType::RAECKE_FRT_FLAT;
    if (s == "raecke_ckr" || s == "ckr" || s == "c")
        return SolverType::RAECKE_CKR_FLAT;
    if (s == "raecke_mst" || s == "random_mst" || s == "rmst" || s == "mst")
        return SolverType::RAECKE_RANDOM_MST_FLAT;
    if (s == "raecke_frt_mendel" || s == "frt_mendel")
        return SolverType::RAECKE_FRT_MENDELSCALING_FLAT;
    if (s == "raecke_ckr_mendel" || s == "ckr_mendel")
        return SolverType::RAECKE_CKR_MENDELSCALING_FLAT;
    if (s == "cohen" || s == "lp" || s == "applegate" || s == "ac" || s == "l")
        return SolverType::LP_APPLEGATE_COHEN;
    if (s == "electrical_parallel" || s == "elec_par" || s == "e_par")
        return SolverType::ELECTRICAL_PARALLEL_BATCHES;
    if (s == "raecke_frt_pointer" || s == "frt_pointer")
        return SolverType::RAECKE_FRT_POINTER;
    if (s == "raecke_ckr_pointer" || s == "ckr_pointer")
        return SolverType::RAECKE_CKR_POINTER;
    if (s == "raecke_mst_pointer" || s == "random_mst_pointer" || s == "rmst_pointer" || s == "mst_pointer")
        return SolverType::RAECKE_RANDOM_MST_POINTER;
    if (s == "raecke_frt_mendel_pointer" || s == "frt_mendel_pointer")
        return SolverType::RAECKE_FRT_MENDELSCALING_POINTER;
    if (s == "raecke_ckr_mendel_pointer" || s == "ckr_mendel_pointer")
        return SolverType::RAECKE_CKR_MENDELSCALING_POINTER;


    if (s == "0") return SolverType::ELECTRICAL;
    if (s == "1") return SolverType::RAECKE_FRT_FLAT;
    if (s == "2") return SolverType::RAECKE_CKR_FLAT;
    if (s == "3") return SolverType::RAECKE_RANDOM_MST_FLAT;
    if (s == "4") return SolverType::LP_APPLEGATE_COHEN;
    if (s == "5") return SolverType::ELECTRICAL_PARALLEL_BATCHES;
    if (s == "6") return SolverType::RAECKE_FRT_MENDELSCALING_FLAT;
    if (s == "7") return SolverType::RAECKE_CKR_MENDELSCALING_FLAT;
    if (s == "8") return SolverType::RAECKE_FRT_POINTER;
    if (s == "9") return SolverType::RAECKE_CKR_POINTER;
    if (s == "10") return SolverType::RAECKE_RANDOM_MST_POINTER;
    if (s == "11") return SolverType::RAECKE_FRT_MENDELSCALING_POINTER;
    if (s == "12") return SolverType::RAECKE_CKR_MENDELSCALING_POINTER;

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

// Parses a comma-separated list of demand model tokens, e.g. "gravity,uniform,gaussian"
inline std::optional<std::vector<DemandModelType>>
parse_demand_model_list(const std::string& s) {
    std::vector<DemandModelType> result;
    size_t start = 0;
    while (true) {
        size_t pos   = s.find(',', start);
        std::string token = (pos == std::string::npos)
                                ? s.substr(start)
                                : s.substr(start, pos - start);
        auto dm = parse_demand_model_token(token);
        if (!dm) return std::nullopt;
        result.push_back(*dm);
        if (pos == std::string::npos) break;
        start = pos + 1;
    }
    if (result.empty()) return std::nullopt;
    return result;
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
       << "  " << "./oblivious_routing " << " <solver> <graph_file> [<demand_models>] [<graph_format>]\n\n"
       << "Solvers (case-insensitive, comma-separated list allowed):\n"
       << "  electrical | ef | e           -> Electrical Flow (naive)\n"
       << "  raecke_frt | frt | f          -> Tree-based (Raecke/FRT)\n"
       << "  raecke_ckr | ckr | c          -> Tree-based (Raecke/CKR)\n"
       << "  raecke_frt_mendel | frt_mendel -> Tree-based (Raecke/FRT) using MendelScaling\n"
       << "  raecke_ckr_mendel | ckr_mendel -> Tree-based (Raecke/CKR) using MendelScaling\n"
       << "  mst | random_mst | raecke_mst | rmst -> Tree-based (Raecke/Random MST)\n"
       << "  cohen | lp | applegate | ac   -> LP (Applegate and Cohen)\n"
       << "  electrical_parallel | elec_par | e_par -> Electrical Flow (parallel)\n"
       << "Numeric shortcuts: 0=electrical, 1=frt, 2=ckr, 3=mst, 4=LP, 5=electrical_parallel, 6=frt_mendel, 7=ckr_mendel\n"
       << "[Optional] demand models (case-insensitive, comma-separated list allowed):\n"
       << "  gravity | gravity_model       -> Gravity Model\n"
       << "  bimodal | bimodal_model       -> Bimodal Model\n"
       << "  gaussian | gaussian_model     -> Gaussian Model\n"
       << "  uniform | uniform_model       -> Uniform Model\n"
       << "  Example: \"gravity,uniform,gaussian\"\n"
       << "If demand models are provided, the oblivious ratio is compared to the offline optimum for each.\n"
       << "[Optional] graph format (case-insensitive):\n"
       << "  csr | adjlist | list | adjmatrix | matrix\n";
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
        return Config{ *solvers_opt, std::string(argv[2]), {}, GraphFormat::CSR };
    }

    if (argc == 4) {
        // argv[3] can be a demand model list OR a graph format — try demand list first.
        auto demand_opt = parse_demand_model_list(argv[3]);
        if (demand_opt) {
            return Config{ *solvers_opt, std::string(argv[2]), *demand_opt, GraphFormat::CSR };
        }
        auto graph_format_opt = parse_graph_format_token(argv[3]);
        if (graph_format_opt) {
            return Config{ *solvers_opt, std::string(argv[2]), {}, *graph_format_opt };
        }
        if (err) *err = "Unknown demand model list / graph format: " + std::string(argv[3]) + "\n" + usage(argv[0]);
        return std::nullopt;
    }

    if (argc == 5) {
        auto demand_opt = parse_demand_model_list(argv[3]);
        if (!demand_opt) {
            if (err) *err = "Unknown demand model list: " + std::string(argv[3]) + "\n" + usage(argv[0]);
            return std::nullopt;
        }
        auto graph_format_opt = parse_graph_format_token(argv[4]);
        if (!graph_format_opt) {
            if (err) *err = "Unknown graph format: " + std::string(argv[4]) + "\n" + usage(argv[0]);
            return std::nullopt;
        }
        return Config{ *solvers_opt, std::string(argv[2]), *demand_opt, *graph_format_opt };
    }

    if (err) *err = usage(argv[0]);
    return std::nullopt;
}


inline std::string demandModelName(DemandModelType type) {
    switch (type) {
        case DemandModelType::GRAVITY:  return "gravity";
        case DemandModelType::BIMODAL:  return "bimodal";
        case DemandModelType::GAUSSIAN: return "gaussian";
        case DemandModelType::UNIFORM:  return "uniform";
        default:                        return "<unknown>";
    }
}

inline std::unique_ptr<DemandModel> makeDemandModel(DemandModelType type) {
    switch (type) {
        case DemandModelType::GRAVITY:  return std::make_unique<GravityModel>();
        case DemandModelType::BIMODAL:  return std::make_unique<BimodalModel>();
        case DemandModelType::GAUSSIAN: return std::make_unique<GaussianModel>();
        case DemandModelType::UNIFORM:  return std::make_unique<UniformModel>();
        default:
            throw std::runtime_error("Unknown demand model type.");
    }
}

// Calls callback(model_name, demand_map) once per demand model in cfg->demand_models.
inline void HandleDemandModels(const std::optional<Config>& cfg,
                               IGraph& g,
                               std::function<void(const std::string&, const DemandMap&)> callback)
{
    if (!cfg || cfg->demand_models.empty()) return;

    std::vector<std::pair<int,int>> demands;
    demands.reserve(static_cast<size_t>(g.getNumNodes()) * (g.getNumNodes() - 1));
    for (int v : g.getVertices())
        for (int u : g.getVertices())
            if (v != u) demands.push_back({v, u});

    for (DemandModelType type : cfg->demand_models) {
        auto model      = makeDemandModel(type);
        DemandMap dmap  = model->generate(g, demands);
        callback(demandModelName(type), dmap);
    }
}

// Legacy single-map helper — kept for call sites that only need one map.
// Uses the first demand model in the list.
inline void HandleDemandModel(int /*argc*/,
                              char** /*argv*/,
                              const std::optional<Config>& cfg,
                              IGraph& g,
                              DemandMap& demand_map)
{
    if (!cfg || cfg->demand_models.empty()) return;

    std::vector<std::pair<int,int>> demands;
    demands.reserve(static_cast<size_t>(g.getNumNodes()) * (g.getNumNodes() - 1));
    for (int v : g.getVertices())
        for (int u : g.getVertices())
            if (v != u) demands.push_back({v, u});

    auto model = makeDemandModel(cfg->demand_models.front());
    demand_map = model->generate(g, demands);
}

inline void printStatsForDemandModel(const std::string& model_name,
                                     std::pair<double, double> result) {
    if (result.first > 0.0 && result.second > 0.0) {
        std::cout << "Ratio off the optimal offline solution ["
                  << model_name << "] demand model: "
                  << (result.second / result.first) * 100.0 << "% "
                  << "(" << result.first << " / " << result.second << ")\n";
    } else {
        std::cout << "Invalid congestion values for demand model [" << model_name << "].\n";
    }
}

inline double computeRoutingSchemeCongestion(IGraph& _g,
                                             const std::unique_ptr<RoutingScheme>& routing_scheme,
                                             const DemandMap& demand_map) {
    std::vector<double> congestion_per_edge(_g.getNumEdges(), 0.0);
    routing_scheme->routeDemands(congestion_per_edge, demand_map);
    double max_cong = routing_scheme->getMaxCongestion(congestion_per_edge);
    for (const auto& cong : congestion_per_edge)
        if (cong > max_cong) max_cong = cong;
    return max_cong;
}

inline double computeOfflineOptimalCongestion(IGraph& _g, const DemandMap& demand_map) {
    CMMF_Solver mccf(_g);
    mccf.init(_g);
    mccf.AddDemandMap(demand_map);
    auto offline_scheme = mccf.solve();
    return mccf.getCongestion();
}


#endif //OBLIVIOUSROUTING_PARSE_PARAMETER_H
