#ifndef OBLIVIOUSROUTING_ELECTRICAL_FLOW_PAR_BATCHES_H
#define OBLIVIOUSROUTING_ELECTRICAL_FLOW_PAR_BATCHES_H

#include "electrical/electrical.h"

/*
* This is parallel version of the electrical MWU algorithm.
* The parallelism is realized by partitioning the nodes into subsets of ndodes and solving the electrical flow problems for each subset in parallel using a pool of AMG solvers.
*/

class ParallelElectricalMWU : public ElectricalMWU {
public:

    std::vector<std::unique_ptr<LaplacianSolver>> amg_pool;
    int p_threads = 1;

    ParallelElectricalMWU(IGraph& g, int root, bool debug = false):ElectricalMWU(g, root,debug){}
    void run(LinearRoutingTable &table) override;
    void updateEdgeDistances(const std::vector<double> &load) override;
    void getApproxLoad(std::vector<double>& load) override;

    void init( bool debug = false, boost::property_tree::ptree _params = boost::property_tree::ptree()) override;


};


#endif //OBLIVIOUSROUTING_ELECTRICAL_FLOW_PAR_BATCHES_H