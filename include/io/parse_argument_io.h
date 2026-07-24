//
// Created by Mert Biyikli on 25.03.26.
//

#ifndef OBLIVIOUSROUTING_PARSE_ARGURMENT_IO_H
#define OBLIVIOUSROUTING_PARSE_ARGURMENT_IO_H

#include "solver_io.h"
#include "graph_io.h"
#include "demand_io.h"
#include "core/errors.h"
#include "core/config.h"
#include "routing/routing_result.h"
#include <optional>
#include <map>

inline std::string to_lower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c){ return std::tolower(c); });
    return s;
}

// Forward declarations
inline Result<SolverType> parse_solver_token(std::string s);
inline std::optional<DemandModelType> parse_demand_model_token(std::string s);

inline Result<SolverType> parse_solver_token(std::string s) {
    auto it = SOLVER_MAP.find(to_lower(s));
    if  (it != SOLVER_MAP.end()) {
        return (it->second);
    }else {
        return makeErrorMessage(ErrorCode::InvalidSolver, "Solver was not found.");
    }
}

inline std::optional<DemandModelType> parse_demand_model_token(std::string s) {
    auto it = DEMAND_MAP.find(to_lower(std::move(s)));
    if (it != DEMAND_MAP.end())
        return std::optional(it->second);
    else return std::nullopt;
}



// Generic list parser template
template<typename T, typename Parser>
inline Result<std::vector<T>> parse_list(const std::string& s, Parser parser) {
    std::vector<T> result;
    size_t start = 0;
    while (true) {
        size_t pos = s.find(',', start);
        std::string token = (pos == std::string::npos) ? s.substr(start) : s.substr(start, pos - start);
        auto val = parser(token);
        if (!val) return makeErrorMessage(ErrorCode::InvalidArgument, "Input is null.");
        result.push_back(*val);
        if (pos == std::string::npos) break;
        start = pos + 1;
    }
    if  (result.empty()) {
        return makeErrorMessage(ErrorCode::InvalidArgument, "Returned empty input.");
    }
    else {
         return (result);
    }
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
    } catch (const std::exception& /*e*/) { return std::nullopt; } catch (...) { return std::nullopt; }
}

inline std::optional<std::string> parse_output(std::string s) {
    std::string file_name = to_lower(std::move(s));
    if (file_name == "") return std::nullopt;

    auto punc = std::find(file_name.begin(), file_name.end(), '.');
    if (punc != file_name.end()) {
        std::string ext = std::string(punc+1, file_name.end()); // Remove extension if present
        if (ext == "json" || ext == "jason" || ext == "txt") {
            return file_name;
        }
    }

    return std::nullopt;
}


inline std::optional<OutputFormat> parse_output_format(std::string s) {
    s = to_lower(std::move(s));
    if (s == "txt") return OutputFormat::TEXT;
    if (s == "json" || s == "jason") return OutputFormat::JSON;
    if (s == "cout") return OutputFormat::COUT;
    return std::nullopt;
}

inline std::string input_usage() {
    return "Usage:\n <solver> <graph_file> [<demand_models>] [<graph_format>] [<num_threads>] [<output_format>]";
}


// Returns Config on success; prints an error to `err` string on failure.
inline Result<Config> parse_parameter(int argc, char** argv) {
    bool evaluate_demand = false;
    if (argc < 3) {
        return makeErrorMessage(ErrorCode::InputNotFound, ("Too few arguments: " + input_usage()));
    }

    auto solvers_opt = parse_solver_list(argv[1]);
    if (!solvers_opt) {
        return makeErrorMessage(ErrorCode::InputNotFound, "Unknown solver list: " + std::string(argv[1]));
    }

    std::vector<DemandModelType> demands;
    GraphFormat fmt = GraphFormat::CSR;
    OutputFormat out_fmt = OutputFormat::COUT;
    int threads = 1;
    std::string output;

    // Parse optional arguments
    for (int i = 3; i < argc; ++i) {
        std::string arg = argv[i];

        // Try demand models first
        auto d = parse_demand_model_list(arg);
        if (d) { demands = *d; evaluate_demand = true; continue; }

        // Try graph format
        auto g = parse_graph_format_token(arg);
        if (g) { fmt = *g; continue; }

        // Try num threads
        auto t = parse_num_threads(arg);
        if (t) { threads = *t; continue; }


        // Try output
        auto o = parse_output(arg);
        if (o) {output = *o; continue;}

        // Try output format
        auto output_format = parse_output_format(arg);
        if (output_format) { out_fmt = *output_format; continue; }

        return makeErrorMessage(ErrorCode::InvalidArgument, "Unknown argument: " + arg);
    }

    return Config{ *solvers_opt, std::string(argv[2]), evaluate_demand, demands, {}, {},fmt, threads,output, out_fmt};
}




#endif //OBLIVIOUSROUTING_PARSE_ARGURMENT_IO_H