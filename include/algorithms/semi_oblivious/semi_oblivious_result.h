//
// Created by Mert Biyikli on 23.06.26.
//

#ifndef OBLIVIOUSROUTING_SEMI_OBLIVIOUS_RESULT_H
#define OBLIVIOUSROUTING_SEMI_OBLIVIOUS_RESULT_H

struct SemiObliviousRoutingResult {
    std::unique_ptr<RoutingScheme> scheme;
    DemandModelType demand_type{};
    double congestion = -1.0;
    double runtime_microseconds = -1.0;

    std::size_t candidate_paths = 0;
    double average_paths_per_pair = 0.0;
    std::string path_selection_strategy;
};


inline void printSemiObliviousResult(
    const SemiObliviousRoutingResult& r
) {
    std::cout << "Routing base: " << r.path_selection_strategy << std::endl;
    std::cout << "Demand [" << demandModelName(r.demand_type) << "]\n";
    std::cout << "  Congestion: " << r.congestion << '\n';
    std::cout << "  Runtime: " << r.runtime_microseconds << " us\n";
    std::cout << "  Candidate paths: " << r.candidate_paths << '\n';
    std::cout << "  Avg paths/pair: " << r.average_paths_per_pair << '\n';
}

#endif //OBLIVIOUSROUTING_SEMI_OBLIVIOUS_RESULT_H