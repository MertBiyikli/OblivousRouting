//
// Created by Mert Biyikli on 23.10.25.
//

#ifndef OBLIVIOUSROUTING_CKR_PARTITION_H
#define OBLIVIOUSROUTING_CKR_PARTITION_H


#include "../../graph.h"

class CKRPartion {
    Graph m_graph;
public:
    void init(Graph &g, bool debug=false);
    std::vector<int> computePartition(const std::vector<int>& _X, const double& delta);
};

#endif //OBLIVIOUSROUTING_CKR_PARTITION_H