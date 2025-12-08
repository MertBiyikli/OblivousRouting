//
// Created by Mert Biyikli on 25.11.25.
//

#ifndef OBLIVIOUSROUTING_IGRAPH_H
#define OBLIVIOUSROUTING_IGRAPH_H
#include <vector>
#include <fstream>
#include <sstream>
#include <numeric>


#define INVALID_EDGE_ID -1

using Edge = std::pair<int, int>; // (u, v)


// Optionally, define a custom comparator to enforce u < v ordering
struct EdgeCompare {
    bool operator()(const Edge& a, const Edge& b) const {
        return (a.first == b.first) ? (a.second < b.second) : (a.first < b.first);
    }
};




class IGraph {
public:
    int n = 0, m = 0; // Number of nodes and edges
    bool is_processed = false;
    std::vector<int> vertices; // List of vertices

    IGraph() = default;

    explicit IGraph(int n) : n(n), m(0) {
        vertices.resize(n);
        std::iota(vertices.begin(), vertices.end(), 0);
    }

    IGraph(const IGraph&) = default;
    IGraph(IGraph&&) noexcept = default;
    IGraph& operator=(const IGraph&) = default;
    IGraph& operator=(IGraph&&) noexcept = default;

    virtual ~IGraph() = default;

    virtual void finalize() = 0;

    struct NeighborRangeBase {
        virtual ~NeighborRangeBase() = default;
        virtual const int* beginPtr() const = 0;
        virtual const int* endPtr() const = 0;
    };

    struct NeighborRange {
        const int* b;
        const int* e;

        const int* begin() const { return b; }
        const int* end()   const { return e; }

        const int size() const { return static_cast<int>(e - b); }
        int operator[](size_t i) const {
            return b[i];
        }
    };

    virtual NeighborRange neighbors(int node) const = 0;

    // node index based access
    virtual double getEdgeCapacity(int u, int v) const = 0;
    virtual double getEdgeDistance(int u, int v) const = 0;
    virtual bool updateEdgeDistance(int u, int v, double distance) = 0;
    virtual int getEdgeId(int u, int v) const = 0;

    // edge index based access
    virtual double getEdgeCapacity(int e) const = 0;
    virtual double getEdgeDistance(int e) const = 0;
    virtual bool updateEdgeDistance(int e, double dist) = 0;
    virtual std::pair<int, int> edgeEndpoints(int e) const = 0;
    virtual const int getAntiEdge(int e) const = 0;

    virtual void addEdge(int u, int v, double capacity, double distance = 1) = 0;
    virtual std::vector<int> getShortestPath(int source, int target) const = 0;
    virtual const double GetDiameter() const = 0;

/*
 *The following methods are for the parameterized shortest path computations.
 *These methods can be overridden in derived classes if needed.
 */
    template<typename T>
    T getDistanceVector() {
        T dist;
        // Special-case adjacency-list style distances
        if constexpr (std::is_same_v<T, std::vector<std::vector<double>>>) {
            dist.resize(n);
            for (int u = 0; u < n; ++u) {
                dist[u].assign(neighbors(u).size(), 1.0); // default distance 1.0
                for (auto& v : neighbors(u)) {
                    dist[u].emplace_back(getEdgeDistance(u, v));
                }
            }
            return dist;
        }
        // CSR style
        else if constexpr (std::is_same_v<T, std::vector<double>>) {
            dist.resize(m);
            for (int e = 0; e < m; ++e) {
                dist[e] = getEdgeDistance(e);
            }
            return dist;
        }
        // Fallback: default-construct
        else {
            return T{};
        }
    }

    virtual std::vector<int>
    getShortestPath(int s, int t, const std::vector<double>& dist) const = 0;



    virtual int getNumNodes() const {
        return n;
    }

    virtual int getNumEdges() const {
        return m;
    }

    virtual const std::vector<int>& getVertices() const {
        return vertices;
    }

    virtual void resetEdgeWeights() {
        for (int v = 0; v < n; ++v) {
            for (int i = 0; i < neighbors(v).size(); ++i) {
                updateEdgeDistance(v, neighbors(v)[i], 1.0);
            }
        }
    }

    // Note that this function does change the members from the base class
    virtual void InitializeMemberByParser(int maxNodeIdSeen) = 0;

};


inline void readLGFFile(IGraph& g, const std::string& filename, bool withDistances) {
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
        g.InitializeMemberByParser(maxNodeIdSeen);


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
            g.addEdge((u), (v), capacityValue);
        }
    }



#endif //OBLIVIOUSROUTING_IGRAPH_H
