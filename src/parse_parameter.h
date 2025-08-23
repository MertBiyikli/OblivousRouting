//
// Created by Mert Biyikli on 23.08.25.
//

#ifndef OBLIVIOUSROUTING_PARSE_PARAMETER_H
#define OBLIVIOUSROUTING_PARSE_PARAMETER_H
#include <string>

enum class SolverType {
    ELECTRICAL_NAIVE,
    RAECKE_FRT,            // tree-based
    LP_APPLEGATE_COHEN
};

struct Config {
    SolverType  solver;
    std::string filename;
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

    // Optional: numeric shortcuts (documented in usage)
    if (s == "0") return SolverType::ELECTRICAL_NAIVE;
    if (s == "1") return SolverType::RAECKE_FRT;
    if (s == "2") return SolverType::LP_APPLEGATE_COHEN;

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
       << "Numeric shortcuts: 0=electric, 1=tree, 2=cohen\n";
    return os.str();
}

// Returns Config on success; prints an error to `err` string on failure.
inline std::optional<Config> parse_parameter(int argc, char** argv, std::string* err) {
    if (argc != 3) {
        if (err) *err = usage(argv[0]);
        return std::nullopt;
    }
    auto solver_opt = parse_solver_token(argv[1]);
    if (!solver_opt) {
        if (err) *err = "Unknown solver: " + std::string(argv[1]) + "\n" + usage(argv[0]);
        return std::nullopt;
    }
    Config cfg{*solver_opt, std::string(argv[2])};
    return cfg;
}
#endif //OBLIVIOUSROUTING_PARSE_PARAMETER_H