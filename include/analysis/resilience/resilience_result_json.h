#ifndef E_ROUTING_RESILIENCE_RESULT_JSON_H
#define E_ROUTING_RESILIENCE_RESULT_JSON_H

#include "analysis/resilience/resilience_analyzer.h"
#include "analysis/scenario_result_json.h" // Reuse numeric/network JSON serialization.
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <locale>
#include <ostream>
#include <stdexcept>
#include <string>
#include <system_error>



class ResilienceResultJsonWriter {
public:
    static void writeTo(std::ostream& os, const N1Result& result) {

        const auto flags = os.flags();
        const auto precision = os.precision();
        const auto locale = os.getloc();
        os.imbue(std::locale::classic());
        os << std::setprecision(std::numeric_limits<double>::max_digits10);
        try {
            std::size_t successful = 0;
            std::size_t disconnected = 0;
            std::size_t execution_failed = 0;
            for (const auto& event : result.outcomes) {
                switch (event.status) {
                    case FailureStatus::SUCCESS:
                        if (!event.reoptimized || !event.comparison)
                            throw std::runtime_error("Successful N-1 event missing measurements");
                        ++successful;
                        break;
                    case FailureStatus::DISCONNECTED:
                        if (event.reoptimized || event.comparison)
                            throw std::runtime_error("Disconnected N-1 event contains measurements");
                        ++disconnected;
                        break;
                    case FailureStatus::EXECUTION_FAILED:
                        if (event.reoptimized || event.comparison)
                            throw std::runtime_error("Failed N-1 event contains measurements");
                        ++execution_failed;
                        break;
                }
            }
            os << "{\"schema_version\":\"1.0\",\"analysis_kind\":\"n1_resilience\",\"graph_fingerprint\":";
            quoted(os, std::to_string(result.graph_fingerprint));
            os << ",\"solver\":";
            quoted(os, getSolverTypeName(result.solver));
            os << ",\"baseline\":";
            network(os, std::optional<NetworkData>{result.baseline});
            os << ",\"summary\":{\"tested_links\":" << result.outcomes.size()
               << ",\"successful\":" << successful
               << ",\"disconnected\":" << disconnected
               << ",\"execution_failed\":" << execution_failed << "}";
            os << ",\"outcomes\":[";
            for (std::size_t i = 0; i < result.outcomes.size(); ++i) {
                if (i) os << ',';
                const auto& event = result.outcomes[i];
                os << "{\"edge_id\":" << event.edge_id
                   << ",\"source\":" << event.source
                   << ",\"target\":" << event.target
                   << ",\"status\":";
                switch (event.status) {
                    case FailureStatus::SUCCESS: quoted(os, "success"); break;
                    case FailureStatus::DISCONNECTED: quoted(os, "disconnected"); break;
                    case FailureStatus::EXECUTION_FAILED: quoted(os, "execution_failed"); break;
                }
                os << ",\"reoptimized\":";
                network(os, event.reoptimized);
                os << ",\"comparison\":";
                comparison(os, event.comparison);
                os << '}';
            }
            os << "]}";
        } catch (...) {
            os.flags(flags);
            os.precision(precision);
            os.imbue(locale);
            throw;
        }
        os.flags(flags);
        os.precision(precision);
        os.imbue(locale);
    }

    static Result<void> write(const N1Result& result, const std::filesystem::path& path) {
        try {
            if (!path.parent_path().empty()) std::filesystem::create_directories(path.parent_path());
            std::ofstream out(path, std::ios::trunc);
            if (!out) return makeErrorMessage(ErrorCode::RuntimeError, "Cannot open N-1 output: " + path.string());
            writeTo(out, result);
            out << '\n';
            out.flush();
            if (!out) return makeErrorMessage(ErrorCode::RuntimeError, "Cannot write N-1 output: " + path.string());
            return {};
        } catch (const std::exception& e) {
            return makeErrorMessage(ErrorCode::RuntimeError, "N-1 serialization failed: " + std::string(e.what()));
        }
    }
};

#endif // E_ROUTING_RESILIENCE_RESULT_JSON_H
