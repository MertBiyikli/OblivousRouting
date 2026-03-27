//
// Created by Mert Biyikli on 20.03.26.
//

#include <sstream>
#include <fstream>
#include <algorithm>
#include <unordered_set>
#include "../../include/io/graph_io.h"
#include "../../include/utils/hash.h"

void readLGFFile(IGraph& g, const std::string& filename) {
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
            } catch (...) {
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
                } catch (...) {
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
                } catch (...) {
                }
            }

            if (!capacityParsed && costParsed && costValue != 0.0) {
                // use inverse of cost as capacity
                capacityValue = 1.0 / costValue;
            }

            if (!capacityParsed && !costParsed) {
                // both capacity and cost parsing failed, use default capacity = 1.0
                capacityValue = 1.0;
            }

            g.addEdge((u), (v), capacityValue);
        }
    }