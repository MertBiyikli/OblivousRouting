
#ifndef OBLIVIOUSROUTING_CONFIG_H// By hand :)
#define OBLIVIOUSROUTING_CONFIG_H

#include "types.h"
#include "utils/demands.h"

struct Config {
    std::vector<SolverType>      solvers;
    std::string                  filename;
    bool evaluate_demand_models = false;
    std::vector<DemandModelType> demand_models;
    std::map<std::string, double> offline_opt_per_model;
    std::map<std::string, demands> demand_maps;
    GraphFormat                  graph_format;
    int num_threads = 1;
    std::string output_filename;
    OutputFormat                 output_format;
    int seed = 42;
};




#endif //OBLIVIOUSROUTING_CONFIG_H