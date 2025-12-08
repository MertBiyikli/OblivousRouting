//
// Created by Mert Biyikli on 23.10.25.
//

#ifndef OBLIVIOUSROUTING_CKR_PARTITION_H
#define OBLIVIOUSROUTING_CKR_PARTITION_H


#include "../../datastructures/GraphADJ.h"
#include "../../datastructures/GraphCSR.h"


struct CKRLevel {
    double R = 0.0;                      // the random radius used at this level
    std::vector<int> owner;              // owner[v] = center that captured v at this level (or -1)
    std::vector<int> pred;               // pred[v] = predecessor of v towards its owner center at this level (-1 for center)
    std::vector<int> centers;            // the centers chosen at this level (subset of V)
};

class CKRPartition {
    GraphADJ m_graph;
public:
    void init(const GraphADJ &g, bool debug=false);

    std::vector<int> computePartition(const std::vector<int>& _X, const double& delta);
    std::vector<int> computePartition(const std::vector<int>& X, const double& delta, CKRLevel& L);
    std::vector<int> computePartition(const GraphCSR& g, const std::vector<int>& X, const double& delta, CKRLevel& L);
};

#endif //OBLIVIOUSROUTING_CKR_PARTITION_H