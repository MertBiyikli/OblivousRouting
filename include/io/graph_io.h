//
// Created by Mert Biyikli on 20.03.26.
//

#ifndef OBLIVIOUSROUTING_GRAPH_IO_H
#define OBLIVIOUSROUTING_GRAPH_IO_H

#include "../data_structures/graph/Igraph.h"
#include "../data_structures/graph/graph_csr.h"
#include "../data_structures/graph/graph_adj.h"


enum class GraphFormat {
    CSR,
    ADJLIST
};

void readLGFFile(IGraph& g, const std::string& filename);


inline std::unique_ptr<IGraph> makegraph(GraphFormat type) {
    switch (type) {
        case GraphFormat::CSR:
            return std::make_unique<GraphCSR>();

        case GraphFormat::ADJLIST:
            return std::make_unique<GraphADJList>();

        default:
            throw std::runtime_error("Unknown graph format.");
    }
}
#endif //OBLIVIOUSROUTING_GRAPH_IO_H