//
// Created by Mert Biyikli on 05.11.25.
//

#ifndef OBLIVIOUSROUTING_GRAPH_CSR_H
#define OBLIVIOUSROUTING_GRAPH_CSR_H
#include <algorithm>
#include <fstream>
#include <numeric>
#include <stdexcept>
#include <utility>
#include <vector>
#include <iostream>
#include <sstream>
#include <boost/container/flat_set.hpp>

#include <map>

#define INVALID_EDGE_ID -1

using _Edge = std::pair<int, int>; // (u, v)

// Optionally, define a custom comparator to enforce u < v ordering
struct EdgeCompare {
    bool operator()(const _Edge& a, const _Edge& b) const {
        return (a.first == b.first) ? (a.second < b.second) : (a.first < b.first);
    }
};

/*
 * The idea is to have a graph class that uses CSR format internally
 * to store edges, capacities and distances for faster access.
 */
class Graph_csr {
public:
    int n, m; // number of nodes and edges
    std::vector<int> vertices; // List of vertices
    std::vector<std::pair<int, int>> edges;
    std::vector<int> head;
    std::vector<double> capacity;
    std::vector<double> distance;

    mutable std::vector<double> dist_buf;
    mutable std::vector<int> parent_buf;

    boost::container::flat_set<_Edge, EdgeCompare> edgeSet;


    // TODO: somehow you also have to explicitly state the the default constructor should not be used
    Graph_csr(int n) : n(n), m(0) {
        vertices.resize(n);
        std::iota(vertices.begin(), vertices.end(), 0);
    }
    void preprocess() {

        // check consistency
        if (edges.size() != capacity.size() || edges.size() != distance.size())
            throw std::runtime_error("edges, capacity and distance size mismatch");

        // 1️⃣ Duplicate edges symmetrically
        int orig_m = static_cast<int>(edges.size());
        edges.reserve(2 * orig_m);
        capacity.reserve(2 * orig_m);
        distance.reserve(2 * orig_m);

        for (int i = 0; i < orig_m; ++i) {
            auto [u, v] = edges[i];
            if (u == v) continue; // skip self-loops if undesired
            edges.emplace_back(v, u);
            capacity.emplace_back(capacity[i]);
            distance.emplace_back(distance[i]);
        }

        m = static_cast<int>(edges.size());

        // 2️⃣ Sort by source vertex
        std::vector<size_t> idx(m);
        std::iota(idx.begin(), idx.end(), 0);
        std::stable_sort(idx.begin(), idx.end(), [&](size_t a, size_t b) {
            const auto &A = edges[a];
            const auto &B = edges[b];
            if (A.first != B.first) return A.first < B.first;
            return A.second < B.second;
        });

        std::vector<std::pair<int,int>> sorted_edges(m);
        std::vector<double> sorted_cap(m);
        std::vector<double> sorted_dist(m);
        for (size_t i = 0; i < m; ++i) {
            sorted_edges[i] = edges[idx[i]];
            sorted_cap[i] = capacity[idx[i]];
            sorted_dist[i] = distance[idx[i]];
        }
        edges.swap(sorted_edges);
        capacity.swap(sorted_cap);
        distance.swap(sorted_dist);

        // 3️⃣ Build head array (n+1 elements)
        head.assign(n + 1, 0);
        int curr = 0;
        for (int u = 0; u < n; ++u) {
            while (curr < m && edges[curr].first == u)
                ++curr;
            head[u + 1] = curr;
        }
    }

    inline auto neighbors(int u) const {
        struct NeighborRange {
            const Graph_csr &g;
            int u;
            auto begin() const { return g.edges.begin() + g.head[u]; }
            auto end() const { return g.edges.begin() + g.head[u + 1]; }
        };
        return NeighborRange{*this, u};
    }

    int getNumNodes() const {
        return n;
    }

    int getNumEdges() const {
        return m;
    }

    void addEdge(int u, int v, double cap) {

        if (u > v) std::swap(u, v);
        _Edge edge{u, v};

        auto [it, inserted] = edgeSet.insert(edge);
        if (!inserted) {
            return; // Edge already exists
        }

        edges.emplace_back(u, v);
        capacity.push_back(cap);
        distance.push_back(1.0); // default distance
        m++;
    }

    // returns edge id for (u,v), or -1 if not exists
    int getEdgeId(int u, int v) const {
        for (int e = head[u]; e < head[u+1]; ++e) {
            if (edges[e].second == v)
                return e;
        }
        return INVALID_EDGE_ID;
    }

    double getEdgeCapacity(int edge_id) const {
        if (edge_id > m || edge_id < 0) {
            throw std::out_of_range("Edge index out of range");
        }
        if (capacity.size() < edge_id) {
            throw std::runtime_error("Edge index exceeds capacity size");
        }
        return capacity[edge_id];
    }

    double getEdgeCapacity(int u, int v) const {
        if ( u > v) std::swap(u, v);
        int source = u;
        int target = v;

        _Edge k{u, v};
        if (!edgeSet.contains(k))
            throw std::runtime_error("Edge not found");
        int e = getEdgeId(u, v);
        if (e == INVALID_EDGE_ID)
            throw std::runtime_error("Edge not found");
        return capacity[e];
    }

    double getEdgeDistance(int edge_id) const {
        if (edge_id > m || edge_id < 0) {
            throw std::out_of_range("Edge index out of range");
        }
        if (distance.size() < edge_id) {
            throw std::runtime_error("Edge index exceeds distance size");
        }
        return distance[edge_id];
    }

    double getEdgeDistance(int u, int v) const {
        int source = std::min(u, v);
        int target = std::max(u, v);

        int head_idx = head[source];
        int end_idx = (source + 1 < n) ? head[source + 1] : m;
        for (int e = head_idx; e < end_idx; ++e) {
            const auto& [a, b] = edges[e];
            if (a != source) break; // moved past source's edges
            if (b == target) {
                return distance[e];
            }
        }
        throw std::runtime_error("Edge not found");
    }

    void updateDistance(int edge_id, double dist) {
        if (edge_id > m || edge_id < 0) {
            throw std::out_of_range("Edge index out of range");
        }
        if (distance.size() < edge_id) {
            throw std::out_of_range("Edge index exceeds distance size");
        }
        distance[edge_id] = dist;
    }

    void updateEdgeDistance(int u, int v, double dist) {
        int source = std::min(u, v);
        int target = std::max(u, v);

        if (source < 0 || source >= n || target < 0 || target >= n) {
            throw std::out_of_range("Node index out of range");
        }

         _Edge k{u, v};
        if (!edgeSet.contains(k))
            throw std::runtime_error("Edge not found");


        int head_idx = head[source];
        int end_idx = (source + 1 < n) ? head[source + 1] : m;
        for (int e = head_idx; e < end_idx; ++e) {
            const auto& [a, b] = edges[e];
            if (a != source) break; // moved past source's edges
            if (b == target) {
                distance[e] = dist;
                return;
            }
        }
        throw std::runtime_error("Edge not found");
    }

    std::vector<int> getShortestPath(int src, int tgt = -1, bool full_search = false) const {
        using P = std::pair<double, int>; // (distance, node)
        const double INF = std::numeric_limits<double>::infinity();

        // 🔧 Ensure reusable buffers are the right size
        if ((int)dist_buf.size() != n) {
            dist_buf.resize(n);
            parent_buf.resize(n);
        }

        // 🔁 Reset only what’s needed
        std::fill(dist_buf.begin(), dist_buf.end(), INF);
        std::fill(parent_buf.begin(), parent_buf.end(), -1);

        auto &dist = dist_buf;
        auto &parent = parent_buf;

        dist[src] = 0.0;

        // min-heap priority queue
        std::priority_queue<P, std::vector<P>, std::greater<>> pq;
        pq.emplace(0.0, src);

        while (!pq.empty()) {
            auto [du, u] = pq.top();
            pq.pop();

            if (du > dist[u]) continue;     // outdated entry
            if (u == tgt) break;            // early stop when target reached

            // ✅ early stop only if not full search and target reached
            if (!full_search && tgt >= 0 && u == tgt)
                break;

            // iterate over neighbors via CSR
            for (int e = head[u]; e < head[u + 1]; ++e) {
                int v = edges[e].second;
                double w = distance[e];
                double nd = du + w;
                if (nd < dist[v]) {
                    dist[v] = nd;
                    parent[v] = u;
                    pq.emplace(nd, v);
                }
            }
        }
        // If we don't have a specific target, skip reconstruction
        if (tgt < 0 || full_search)
            return {};
        // reconstruct path (if reachable)
        std::vector<int> path;
        if (dist[tgt] == INF) return path; // no path

        for (int v = tgt; v != -1; v = parent[v])
            path.push_back(v);
        std::reverse(path.begin(), path.end());
        return path;
    }

    // ===============================================================
    // 🔹 Extract edge IDs along the path found by the last Dijkstra run
    // ===============================================================
    std::vector<int> getPathEdgesFromParent(int src, int tgt) const {
        std::vector<int> edges_on_path;
        if (tgt < 0 || tgt >= n) return edges_on_path;
        if (parent_buf.empty()) return edges_on_path;

        for (int v = tgt; v != src && v != -1; v = parent_buf[v]) {
            int u = parent_buf[v];
            if (u == -1) break;

            // Find edge ID for (u,v)
            int start = head[u];
            int end = head[u + 1];
            for (int e = start; e < end; ++e) {
                if (edges[e].second == v) {
                    edges_on_path.push_back(e);
                    break;
                }
            }
        }

        std::reverse(edges_on_path.begin(), edges_on_path.end());
        return edges_on_path;
    }

    // ===============================================================
    // 🔹 Compute shortest path and return edge indices along it
    // ===============================================================
    std::vector<int> getPathEdges(int src, int tgt) const {
        getShortestPath(src, tgt);               // fills dist_buf and parent_buf
        return getPathEdgesFromParent(src, tgt); // reconstruct via parent links
    }




    double getShortestDistance(int src, int tgt, bool full_search = false) const {
        const double INF = std::numeric_limits<double>::infinity();

        // run Dijkstra (fills dist_buf)
        getShortestPath(src, tgt, full_search);

        // ensure target index is valid
        if (tgt < 0 || tgt >= n) return INF;

        // return the precomputed distance value
        return dist_buf[tgt];
    }


    double computeExactDiameterUsingDijkstra() const {
        if (n == 0) return 0.0;

        double diameter = 0.0;

        for (int src = 0; src < n; ++src) {
            // We only need the distances, not the path to any particular target.
            // Passing tgt = src just triggers Dijkstra from src (the early-exit won't trigger).
            getShortestPath(src, -1, true);  // full Dijkstra from src

            // reuse the existing dist_buf from the last run
            const auto &dist = dist_buf;

            // find the farthest reachable node
            double local_max = 0.0;
            for (double d : dist)
                if (d < std::numeric_limits<double>::infinity())
                    local_max = std::max(local_max, d);

            diameter = std::max(diameter, local_max);
        }

        return diameter;
    }



    void print() {
        for (int e = 0; e<m; e++) {
            std::cout << "Edge " << e << " =("<< edges[e].first << ", "  << edges[e].second << " cap:" << capacity[e] << " dist: " << distance[e] << ")" << std::endl;
        }

        // print head array
        std::cout << "Head array: ";
        for (int i = 0; i < n; i++) {
            std::cout << head[i] << " ";
        }
    }

    void readLFGFile(const std::string &filename, bool withDistances) {
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
    n = maxNodeIdSeen + 1;
    m = 0;
        vertices.resize(n);
    std::iota(vertices.begin(), vertices.end(), 0); // fill vertices with 0, 1, ..., maxNodeIdSeen

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
        this->addEdge((u), (v), capacityValue);
    }

}
};

#endif //OBLIVIOUSROUTING_GRAPH_CSR_H