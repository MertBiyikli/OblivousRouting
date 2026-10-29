
#ifndef OBLIVIOUSROUTING_CONFIG_H// By hand :)
#define OBLIVIOUSROUTING_CONFIG_H

#include "types.h"
#include "utils/demands.h"


struct _CapacityUpgrade {
    int edge_id = -1;
    double factor = 1.0;
};

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

    bool failure_analysis = false;
    bool failure_recovery = false;


    std::string routing_table_input_file{};

    std::string routing_table_output_directory{};
    /*
     * Empty means that visualization export is disabled.
     */
    std::string visualization_output_directory{};


    // Inside Config:
    double demand_scale = 1.0;
    std::vector<_CapacityUpgrade> capacity_upgrades;
    std::vector<int> failed_links;

};






#endif //OBLIVIOUSROUTING_CONFIG_H