//
// Created by Mert Biyikli on 20.03.26.
//

#ifndef OBLIVIOUSROUTING_GRAPH_IO_H
#define OBLIVIOUSROUTING_GRAPH_IO_H

#include "../data_structures/graph/Igraph.h"
#include "../data_structures/graph/graph_csr.h"
#include "../data_structures/graph/graph_adj.h"
#include "core/errors.h"
#include "core/types.h"
#include "core/config.h"
#include <memory>


class GraphIO {
public:
    static Result<void> readLGFFile(IGraph& g, const std::string& filename) {
        std::ifstream infile(filename);
        if (!infile) {
            return std::unexpected(Error{ErrorCode::FileNotFound, "Could not open LGF file: " + filename});
        }

        std::string line;
        bool inNodesSection = false;
        bool inArcsSection = false;

        int maxNodeIdSeen = -1;
        int costColIdx = -1;
        int capColIdx  = -1;

        // Phase 1: Read lines until we finish @nodes and find @arcs header.
        std::vector<std::string> allLines;
        while (std::getline(infile, line)) {
            // Trim leading/trailing whitespace:
            // (Here we just check if line is empty or starts with a comment character)
            if (line.empty()) {
                continue;
            }
            allLines.push_back(line);
        }
        infile.close();

        try {
            // 1) First pass: scan for "@nodes" → read node lines → track maxNodeIdSeen.
            for (size_t i = 0; i < allLines.size(); ++i) {
                const std::string &raw = allLines[i];
                std::string lower = raw;
                // normalize to lowercase to match "@nodes" or "@arcs"
                std::transform(lower.begin(), lower.end(), lower.begin(),
                               [](unsigned char c) { return std::tolower(c); });

                if (!inNodesSection) {
                    if (lower.rfind("@nodes", 0) == 0) {
                        inNodesSection = true;
                        continue;
                    }
                } else {
                    // We are inside the @nodes block. Stop if we see "@arcs".
                    if (lower.rfind("@edges", 0) == 0 ||
                        lower.rfind("@arcs", 0) == 0) {
                        inNodesSection = false;
                        inArcsSection = true;
                        // The next line after this is expected to be the arcs‐header.
                        continue;
                        }
                    // Otherwise, we expect this line to have: "<label><TAB><node_id>"
                    // We split on tabs:
                    std::istringstream iss(raw);
                    std::vector<std::string> tokens;
                    std::string tok;
                    while (iss >> tok) {
                        tokens.push_back(tok);
                    }
                    if (tokens.size() >= 2) {
                        // node_id is the second column
                        try {
                            int nid = std::stoi(tokens[1]);
                            maxNodeIdSeen = std::max(maxNodeIdSeen, nid);
                        } catch (const std::exception& /*e*/) {
                            // malformed line: skip (e.what() available if logging)
                        }
                    }

                    try {
                        int nid = std::stoi(tokens[0]);
                        maxNodeIdSeen = std::max(maxNodeIdSeen, nid);
                    } catch (const std::exception& /*e*/) {
                        // malformed line: skip (e.what() available if logging)
                    }

                }
            }

            if (maxNodeIdSeen < 0) {
                return makeErrorMessage(ErrorCode::InvalidGraph,"LGF file had no valid @nodes section or no node IDs found");
            }

            // 2) Initialize our graph with (maxNodeIdSeen+1) nodes
            g.initializeMemberByParser(maxNodeIdSeen);


            // 3) Second pass: scan for "@arcs", find header, then parse each arc line:
            inArcsSection = false;
            bool readHeader = false;
            std::unordered_set<std::pair<int, int>, PairHash> existingEdges;
            bool capacityInHeader = false;
            bool costInHeader = false;
            for (size_t i = 0; i < allLines.size(); ++i) {
                const std::string &raw = allLines[i];
                std::string lower = raw;
                std::transform(lower.begin(), lower.end(), lower.begin(),
                               [](unsigned char c) { return std::tolower(c); });

                if (!inArcsSection) {
                    // Look for “@arcs”
                    if (lower.rfind("@edges", 0) == 0
                        || lower.rfind("@arcs", 0) == 0) {
                        inArcsSection = true;
                        }
                    continue;
                }

                // We are in @arcs section. The very first non‐empty line (after "@arcs") is the header.
                if (inArcsSection && !readHeader) {
                    // This is the header line, e.g. "    label    cost    capacity"
                    std::istringstream iss(raw);
                    std::vector<std::string> headerTokens;
                    std::string tok;
                    while (iss >> tok) {
                        if(!tok.empty()) {
                            headerTokens.push_back(tok);
                        }
                    }
                    // Find which index is “cost” and which is “capacity”
                    bool labelInHeader = false;
                    for (int col = 0; col < (int)headerTokens.size(); ++col) {
                        std::string h = headerTokens[col];
                        std::transform(h.begin(), h.end(), h.begin(),
                                       [](unsigned char c) { return std::tolower(c); });
                        if (h == "cost") {
                            costColIdx = col;
                            costInHeader = true;
                        }
                        if (h == "capacity") {
                            capColIdx = col;
                            capacityInHeader = true;
                        }

                        // if label is also in the header, the column indices shift by 3 (u, v, label), otherwise it would be 2 (u, v)
                        if (h == "label") {
                            labelInHeader = true;
                        }
                    }

                    costColIdx += (costInHeader ? 2 : 0);
                    capColIdx +=  (capacityInHeader ? 2 : 0);
                    /*
                    if (labelInHeader) {
                        // label is header, we need to shift further by 1
                        costColIdx += (costInHeader ? 1 : 0);
                        capColIdx +=  (capacityInHeader ? 1 : 0);
                    }*/
                    readHeader = true;
                    continue;
                }

                // Now each subsequent line (after header) is “u<TAB>v<TAB>…”

                // skip blank or comment lines
                if (raw.empty()
                    || raw[0] == '#') {
                    continue;
                }

                // Split on tabs
                std::istringstream iss(raw);
                std::vector<std::string> parts;
                std::string field;

                while (iss >> field) {
                    if(!field.empty()) {
                        parts.push_back(field);
                    }
                }

                if (parts.size() < 2) {
                    // not enough columns ⇒ skip
                    continue;
                }

                // Parse u, v (1‐based IDs in LGF)
                int u, v;
                try {
                    u = std::stoi(parts[0]);  // convert to 0‐based
                    v = std::stoi(parts[1]);
                } catch (const std::exception& /*e*/) {
                    continue; // malformed
                }
                if (u < 0 || v < 0 || u >= g.n || v >= g.n) {
                    // invalid node Id ⇒ skip
                    continue;
                }

                // We only add each undirected edge once (u < v)
                if (u > v) {
                    std::swap(u, v);
                }

                if (existingEdges.contains({u, v})) {
                    // edge already exists ⇒ skip
                    continue;
                }else {
                    existingEdges.insert({u,v});
                }

                /* Note the following:
                *   in many dataset, such as BackBone, the 'cost' is only provided,
                *   and 'capacity' is missing. In such cases, we choose to
                *   use the inverse of cost as capacity.
                *
                *   if both 'cost' and 'capacity' are provided, we use 'capacity' directly.
                *
                *   if none is provided, we set capacity = 1.0 by default.
                 */

                bool capacityParsed = false;
                double capacityValue = 1.0;
                if (capacityInHeader
                    && capColIdx >= 0
                    && capColIdx < (int)parts.size()) {
                    // attempt to parse capacity
                    try {
                        capacityValue = std::stod(parts[capColIdx]);
                        capacityParsed = true;
                    } catch (const std::exception& /*e*/) {
                        // failed to parse capacity, will try cost or default below
                    }
                }


                bool costParsed = false;
                double costValue = 1.0;
                if (costInHeader
                    && costColIdx >= 0
                    && costColIdx < (int)parts.size()) {
                    // attempt to parse capacity
                    try {
                        costValue = std::stod(parts[costColIdx]);
                        costParsed = true;
                    } catch (const std::exception& /*e*/) {
                    }
                }

                if (!capacityParsed && costParsed && costValue != 0.0) {
                    // use inverse of cost as capacity
                    capacityValue =  costValue;
                }

                if (!capacityParsed && !costParsed) {
                    // both capacity and cost parsing failed, use default capacity = 1.0
                    capacityValue = 1.0;
                }

                g.addEdge((u), (v), capacityValue);
            }
            return {};
        } catch (const std::exception& e) {
            return fromStdException(e, ErrorCode::InvalidGraph);
        } catch (...) {
            return fromUnknownException(ErrorCode::InvalidGraph);
        }
    }
};

inline Result<std::unique_ptr<IGraph>> makegraph(GraphFormat type) {
    switch (type) {
        case GraphFormat::CSR:
            return std::make_unique<GraphCSR>();

        case GraphFormat::ADJLIST:
            return std::make_unique<GraphADJList>();

        default:
            return std::unexpected(Error{ErrorCode::InvalidGraph, "Unknown graph format."});
    }
}


inline Result<std::unique_ptr<IGraph>> load_graph(Config& cfg, int argc, char** argv) {
    // Load or create graph
    auto graph = makegraph(cfg.graph_format);
    if ( !graph ){
        return std::unexpected(graph.error());
    }

    if (!cfg.filename.empty()) {
        auto lgf = GraphIO::readLGFFile(*graph.value(), cfg.filename);
        if (!lgf) {
            return std::unexpected(lgf.error());
        }
    }else {
        return std::unexpected(Error{
                ErrorCode::FileNotFound,
                "Parsed filename is empty."
            });
    }

    graph.value()->finalize();
    return std::move(graph.value());
}

inline Result<std::unique_ptr<optimized::Graph<EdgeData>>> load_graph_optimized(Config& cfg, int argc, char** argv) {
    // Load using IGraph first
    auto igraph = load_graph(cfg, argc, argv);
    if (!igraph) {
        return std::unexpected(igraph.error());
    }
    
    // Convert to optimized::Graph
    int n = igraph.value()->getNumNodes();
    std::vector<optimized::Graph<EdgeData>::InputEdge> input_edges;
    
    // Add edges from IGraph
    for (int u = 0; u < n; ++u) {
        for (int v : igraph.value()->neighbors(u)) {
            if (u < v) {
                double weight = igraph.value()->getEdgeDistance(u, v);
                double capacity = igraph.value()->getEdgeCapacity(u, v);
                input_edges.push_back({u, v, EdgeData{capacity, weight}});
            }
        }
    }
    
    auto opt_graph = std::make_unique<optimized::Graph<EdgeData>>(n, input_edges);
    return opt_graph;
}
#endif //OBLIVIOUSROUTING_GRAPH_IO_H