//
// Created by Mert Biyikli on 11.05.25.
//


#include "graph.h"
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <sstream>
#include <queue>
#include <limits>


void Graph::readGraph(const std::string &filename) {
    std::ifstream infile(filename);
    std::cout << "Reading graph from file: " << filename << std::endl;
    if (!infile) {
        throw std::runtime_error("Could not open file: " + filename);
    }

    std::string line;
    int n = 0, m = 0;

    while (std::getline(infile, line)) {
        if (line.empty() || line[0] == 'c') continue;

        std::istringstream iss(line);
        if (line[0] == 'p') {
            std::string tmp;
            iss >> tmp >> tmp >> n >> m;
            break;
        }
    }

    if (n == 0) {
        throw std::runtime_error("Invalid or missing problem line in DIMACS file");
    }


    this->m_adj.resize(n);
    this->m_adj_capacities.resize(n);
    this->m_adj_distances.resize(n);
    this->m_vertices.resize(n);

    // Fill vertices with 0, 1, ..., maxNodeIdSeen
    std::iota(this->m_vertices.begin(), this->m_vertices.end(), 0);

    do {
        if (line.empty() || line[0] == 'c') continue;

        std::istringstream iss(line);
        if (line[0] == 'a') {
            char dummy;
            int u, v;
            double cap = 1.0;
            iss >> dummy >> u >> v;
            if (!(iss >> cap)) cap = 1.0;
            this->addEdge(u - 1, v - 1, cap);  // convert to 0-based
        }
    } while (std::getline(infile, line));
}

void Graph::readLFGFile(const std::string &filename, bool withDistances) {
    std::ifstream infile(filename);
    if (!infile) {
        throw std::runtime_error("Could not open LGF file: " + filename);
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
                } catch (...) {
                    // malformed line: skip
                }
            }

            try {
                int nid = std::stoi(tokens[0]);
                maxNodeIdSeen = std::max(maxNodeIdSeen, nid);
            } catch (...) {
                // malformed line: skip
            }

        }
    }

    if (maxNodeIdSeen < 0) {
        throw std::runtime_error("LGF file had no valid @nodes section or no node IDs found");
    }

    // 2) Initialize our graph with (maxNodeIdSeen+1) nodes
    int n = maxNodeIdSeen + 1;
    int m = 0;
    m_adj.clear();
    m_adj.resize(n);
    m_adj_distances.resize(n);

    m_adj_capacities.resize(n);
    m_vertices.resize(n);
    std::iota(m_vertices.begin(), m_vertices.end(), 0); // fill vertices with 0, 1, ..., maxNodeIdSeen

    // 3) Second pass: scan for "@arcs", find header, then parse each arc line:
    inArcsSection = false;
    bool readHeader = false;
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
            for (int col = 0; col < (int)headerTokens.size(); ++col) {
                std::string h = headerTokens[col];
                std::transform(h.begin(), h.end(), h.begin(),
                               [](unsigned char c) { return std::tolower(c); });
                if (h == "cost") {
                    costColIdx = col;
                }
                if (h == "capacity") {
                    capColIdx = col;
                }
            }
            costColIdx += 2; // TODO: this is hardcoded for now, should be dynamic
            readHeader = true;
            continue;
        }

        // Now each subsequent line (after header) is “u<TAB>v<TAB>…”
        if (raw.empty() || raw[0] == '#') {
            // skip blank or comment lines
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
            u = std::min(std::stoi(parts[0]), std::stoi(parts[1]));  // convert to 0‐based
            v = std::max(std::stoi(parts[0]), std::stoi(parts[1]));
        } catch (...) {
            continue; // malformed
        }
        if (u < 0 || v < 0 || u >= n || v >= n) {
            // invalid node Id ⇒ skip
            continue;
        }
        // We only add each undirected edge once (u < v)
        if (u > v) {
            continue;
        }


        // Determine capacity:
        double capacityValue = 1.0;
        bool capacityParsed = false;
        if (costColIdx >= 0 && costColIdx < (int)parts.size()) {
            // attempt to parse capacity
            try {
                capacityValue = std::stod(parts[costColIdx]);
                capacityParsed = true;
            } catch (...) {
                capacityParsed = false;
            }
        }

        // TODO: for now we only consider the undirected case...
        if (false && !capacityParsed && withDistances && costColIdx >= 0 && costColIdx < (int)parts.size()) {
            // parse cost, then invert it
            try {
                double costVal = std::stod(parts[costColIdx]);
                if (costVal != 0.0) {
                    capacityValue = costVal;
                } else {
                    capacityValue = 1.0;
                }
                capacityParsed = true;
            } catch (...) {
                capacityParsed = false;
            }
        }
        // if still not parsed, we leave capacityValue = 1.0 by default

        if (capacityValue == 0) {
            capacityValue = 1;
        }

        // Finally, add the undirected edge (Graph::addEdge adds both directions)
        // add random edge capacities
        // capacityValue = rand()%this->getNumNodes()+10; // random factor in [0.5, 1.5]
        this->addEdge(u, v, capacityValue);
    }

    // 4) (Optional) If you have a PruneGraph class, prune degree‐one vertices
    //     new PruneGraph().pruneDegreeOne(*this);
    //
    // If you do not have PruneGraph or you don’t want to prune, simply comment out the line above.
}

