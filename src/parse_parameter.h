//
// Created by Mert Biyikli on 23.08.25.
//

#ifndef OBLIVIOUSROUTING_PARSE_PARAMETER_H
#define OBLIVIOUSROUTING_PARSE_PARAMETER_H
#include <string>

enum class SolverType {
    ELECTRICAL_NAIVE,
    RAECKE_FRT,            // tree-based
    LP_APPLEGATE_COHEN,
    RAECKE_RANDOM_MST
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

    // Optional: numeric shortcuts (documented in usage)
    if (s == "0") return SolverType::ELECTRICAL_NAIVE;
    if (s == "1") return SolverType::RAECKE_FRT;
    if (s == "2") return SolverType::LP_APPLEGATE_COHEN;
    if (s == "3") return SolverType::RAECKE_RANDOM_MST;

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
       << "  cohen | lp | applegate | ac   -> LP (Applegate–Cohen)\n"
       << "  mst | random_mst | raecke_mst | rmst -> LP (Raecke/Random MST)\n\n"
       << "Numeric shortcuts: 0=electric, 1=tree, 2=cohen\n"
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
#endif //OBLIVIOUSROUTING_PARSE_PARAMETER_H