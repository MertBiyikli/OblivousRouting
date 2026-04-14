//
// Created by Mert Biyikli on 25.03.26.
//

#ifndef OBLIVIOUSROUTING_PARSE_ARGURMENT_IO_H
#define OBLIVIOUSROUTING_PARSE_ARGURMENT_IO_H

#include "solver_io.h"
#include "graph_io.h"
#include "demand_io.h"

#include <optional>
#include <map>

struct Config {
    std::vector<SolverType>      solvers;
    std::string                  filename;
    std::vector<DemandModelType> demand_models;
    GraphFormat                  graph_format;
    int num_threads = 1;
};

inline std::string to_lower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c){ return std::tolower(c); });
    return s;
}

// Forward declarations
inline std::optional<SolverType> parse_solver_token(std::string s);
inline std::optional<DemandModelType> parse_demand_model_token(std::string s);

inline std::optional<SolverType> parse_solver_token(std::string s) {
    auto it = SOLVER_MAP.find(to_lower(std::move(s)));
    return (it != SOLVER_MAP.end()) ? std::optional(it->second) : std::nullopt;
}

inline std::optional<DemandModelType> parse_demand_model_token(std::string s) {
    auto it = DEMAND_MAP.find(to_lower(std::move(s)));
    return (it != DEMAND_MAP.end()) ? std::optional(it->second) : std::nullopt;
}



// Generic list parser template
template<typename T, typename Parser>
inline std::optional<std::vector<T>> parse_list(const std::string& s, Parser parser) {
    std::vector<T> result;
    size_t start = 0;
    while (true) {
        size_t pos = s.find(',', start);
        std::string token = (pos == std::string::npos) ? s.substr(start) : s.substr(start, pos - start);
        auto val = parser(token);
        if (!val) return std::nullopt;
        result.push_back(*val);
        if (pos == std::string::npos) break;
        start = pos + 1;
    }
    return result.empty() ? std::nullopt : std::optional(result);
}

// Unified list parsers using generic template
inline auto parse_solver_list(const std::string& s) {
    return parse_list<SolverType>(s, parse_solver_token);
}

inline auto parse_demand_model_list(const std::string& s) {
    return parse_list<DemandModelType>(s, parse_demand_model_token);
}

inline std::optional<GraphFormat> parse_graph_format_token(std::string s) {
    s = to_lower(std::move(s));
    if (s == "csr") return GraphFormat::CSR;
    if (s == "adjlist" || s == "list") return GraphFormat::ADJLIST;
    return std::nullopt;
}

inline std::optional<int> parse_num_threads(std::string s) {
    try {
        int num = std::stoi(s);
        return (num > 0) ? std::optional(num) : std::nullopt;
    } catch (...) { return std::nullopt; }
}

// Generate all-pairs demand list
inline std::vector<std::pair<int,int>> generateAllDemandPairs(IGraph& g) {
    std::vector<std::pair<int,int>> result;
    result.reserve(static_cast<size_t>(g.getNumNodes()) * (g.getNumNodes() - 1));
    for (int v : g.getVertices())
        for (int u : g.getVertices())
            if (v != u) result.push_back({v, u});
    return result;
}

// Demand model handling
inline void HandleDemandModels(const std::optional<Config>& cfg, IGraph& g,
                               std::function<void(const std::string&, const demands&)> callback) {
    if (!cfg || cfg->demand_models.empty()) return;
    auto pairs = generateAllDemandPairs(g);
    for (DemandModelType type : cfg->demand_models) {
        auto model = makeDemandModel(type);
        demands dmap = model->generate(g, pairs);
        callback(demandModelName(type), dmap);
    }
}

inline demands GetSingleDemandModel(const std::optional<Config>& cfg, IGraph& g) {
    if (!cfg || cfg->demand_models.empty()) return demands{};
    auto pairs = generateAllDemandPairs(g);
    auto model = makeDemandModel(cfg->demand_models.front());
    return model->generate(g, pairs);
}


inline std::string usage(const char* prog) {
    std::ostringstream os;
    os << "Usage:\n  " << prog << " <solver> <graph_file> [<demand_models>] [<graph_format>] [<num_threads>] [<cycle_strategy>]\n\n"
       << "Solvers (case-insensitive, comma-separated):\n"
       << "  electrical | ef | e | 0                          -> Electrical Flow\n"
       << "  raecke_frt | frt | f | 1                         -> Tree-based (Raecke/FRT)\n"
       << "  raecke_ckr | ckr | c | 2                         -> Tree-based (Raecke/CKR)\n"
       << "  mst | random_mst | rmst | 3                      -> Tree-based (Random MST)\n"
       << "  cohen | lp | applegate | ac | 4                  -> LP (Applegate-Cohen)\n"
       << "  raecke_frt_mendel | frt_mendel | 6               -> FRT with MendelScaling\n"
       << "  raecke_ckr_mendel | ckr_mendel | 7               -> CKR with MendelScaling\n"
       << "  *_pointer variants (8-12)                        -> Pointer-based HST versions\n"
       << "Demand Models (optional, comma-separated):\n"
       << "  gravity | bimodal | gaussian | uniform\n"
       << "Graph Format (optional): csr | adjlist\n";
    return os.str();
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

    std::vector<DemandModelType> demands;
    GraphFormat fmt = GraphFormat::CSR;
    int threads = 1;

    // Parse optional arguments
    for (int i = 3; i < argc; ++i) {
        std::string arg = argv[i];

        // Try demand models first
        auto d = parse_demand_model_list(arg);
        if (d) { demands = *d; continue; }

        // Try graph format
        auto g = parse_graph_format_token(arg);
        if (g) { fmt = *g; continue; }

        // Try num threads
        auto t = parse_num_threads(arg);
        if (t) { threads = *t; continue; }


        if (err) *err = "Unknown argument: " + arg + "\n" + usage(argv[0]);
        return std::nullopt;
    }

    return Config{ *solvers_opt, std::string(argv[2]), demands, fmt, threads};
}


#endif //OBLIVIOUSROUTING_PARSE_ARGURMENT_IO_H