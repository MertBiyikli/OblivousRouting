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



//////////////////////////////////
// Base class methods
//////////////////////////////////
std::vector<double> Graph::getDistances(int u) const {
    if(u < 0 || u >= m_adj.size()) {
        throw std::out_of_range("Node index out of range");
    }
    return m_distanceMatrix[u];
}


IGraph::NeighborRange Graph::neighbors(int u) const {
    const auto& row = m_adj[u];
    return NeighborRange{
        row.data(),
        row.data() + row.size()
    };
}



void Graph::addEdge(int u , int v, double capacity, double distance) {
    if(u >= m_adj.size() || v >= m_adj.size()) {
        throw std::out_of_range("Node index out of range");
    }
    if(u == v) {
        throw std::invalid_argument("No self-loops allowed (u == v)");
    }

    auto it = std::find(m_adj[u].begin(), m_adj[u].end(), v);
    if(it != m_adj[u].end()) {
        return;
        //throw std::runtime_error("Edge already exists");
    }

    m_adj[u].push_back(v);
    m_adj_capacities[u].push_back(capacity);
    m_adj_distances[u].push_back(distance);

    m_adj[v].push_back(u);
    m_adj_capacities[v].push_back(capacity);
    m_adj_distances[v].push_back(distance);

    // 🔽 NEW: register a single logical edge ID for the undirected edge
    int e = static_cast<int>(m_edgeList.size());
    m_edgeList.emplace_back(u, v);

    m_uvToEdge[{u, v}] = e;
    m_uvToEdge[{v, u}] = e;
    m++;
}

int Graph::getEdgeId(int u, int v) const {
    auto it = m_uvToEdge.find({u, v});
    if (it == m_uvToEdge.end()) {
        return INVALID_EDGE_ID; // or INVALID_EDGE_ID if you use that macro
    }
    return it->second;
}


std::pair<int,int> Graph::edgeEndpoints(int e) const {
    if (e < 0 || e >= static_cast<int>(m_edgeList.size())) {
        throw std::out_of_range("edgeEndpoints: edge id out of range");
    }
    return m_edgeList[e];
}


double Graph::getEdgeCapacity(int u, int v) const {
    double capacity = 0.0;
    if(u >= m_adj.size() || v >= m_adj.size()) {
        throw std::out_of_range("Node index out of range");
    }
    auto it = std::find(m_adj[u].begin(), m_adj[u].end(), v);
    if(it != m_adj[u].end()) {
        int index = std::distance(m_adj[u].begin(), it);
        capacity = m_adj_capacities[u][index];
    } else {
        throw std::runtime_error("Edge not found");
    }
    return capacity;
}


double Graph::getEdgeDistance(int u, int v) const {
    double distance = 0.0;
    if(u >= m_adj.size() || v >= m_adj.size()) {
        throw std::out_of_range("Node index out of range");
    }
    auto it = std::find(m_adj[u].begin(), m_adj[u].end(), v);
    if(it != m_adj[u].end()) {
        int index = std::distance(m_adj[u].begin(), it);
        distance = m_adj_distances[u][index];
    } else {
        throw std::runtime_error("Edge not found");
    }
    return distance;
}

void Graph::updateEdgeDistance(int u, int v, double distance) {
    if(u >= m_adj.size() || v >= m_adj.size()) {
        throw std::out_of_range("Node index out of range");
    }
    auto it = std::find(m_adj[u].begin(), m_adj[u].end(), v);

    if(it != m_adj[u].end()) {
        int index = std::distance(m_adj[u].begin(), it);
        m_adj_distances[u][index] = distance;
        auto anti_edge = std::find(m_adj[v].begin(), m_adj[v].end(), u);
        if(anti_edge != m_adj[v].end()) {
            int anti_index = std::distance(m_adj[v].begin(), anti_edge);
            m_adj_distances[v][anti_index] = distance; // Assuming undirected graph, update both directions
        }
    } else {
        throw std::runtime_error("Edge not found");
    }
}

void Graph::updateEdgeCapacity(int u, int v, double capacity) {
    if(u >= m_adj.size() || v >= m_adj.size()) {
        throw std::out_of_range("Node index out of range");
    }
    auto it = std::find(m_adj[u].begin(), m_adj[u].end(), v);
    if(it != m_adj[u].end()) {
        int index = std::distance(m_adj[u].begin(), it);
        m_adj_capacities[u][index] = capacity;

        // update reverse direction accordingly
        auto anti_edge = std::find(m_adj[v].begin(), m_adj[v].end(), u);
        int anti_index = std::distance(m_adj[v].begin(), anti_edge);
        m_adj_capacities[v][anti_index] = capacity; // Assuming undirected graph, update both directions
    } else {
        throw std::runtime_error("Edge not found");
    }
}

double Graph::getEdgeCapacity(int e) const {
    auto [u, v] = edgeEndpoints(e);
    return getEdgeCapacity(u, v);
}

double Graph::getEdgeDistance(int e) const {
    auto [u, v] = edgeEndpoints(e);
    return getEdgeDistance(u, v);
}

void Graph::updateEdgeDistance(int e, double dist) {
    auto [u, v] = edgeEndpoints(e);
    updateEdgeDistance(u, v, dist);
}


const double& Graph::getShortestDistance(int u, int v) const {
    return m_distanceMatrix[u][v];
}

// maybe think of a way better/faster way of computin the diameter
const double Graph::GetDiameter() const {
    if (m_distanceMatrix.empty()) {
        throw std::runtime_error("Distance matrix is not initialized. Call createDistanceMatrix() first.");
    }

    double diameter = 0.0;
    for (const auto& row : m_distanceMatrix) {
        for (double distance : row) {
            if (distance > diameter && distance < std::numeric_limits<double>::infinity()) {
                diameter = distance;
            }
        }
    }
    return diameter;
}




std::vector<int> Graph::getShortestPath(int u, int v) const {
    if (u < 0 || u >= m_adj.size() || v < 0 || v >= m_adj.size()) {
        throw std::out_of_range("Node index out of range");
    }

    const int n = m_adj.size();
    std::vector<double> dist(n, std::numeric_limits<double>::infinity());
    std::vector<int> prev(n, -1);
    std::vector<bool> visited(n, false);

    dist[u] = 0.0;
    using Node = std::pair<double, int>; // (distance, node)
    std::priority_queue<Node, std::vector<Node>, std::greater<>> pq;
    pq.emplace(0.0, u);

    while (!pq.empty()) {
        auto [d, node] = pq.top();
        pq.pop();

        if (visited[node]) continue;
        visited[node] = true;

        if (node == v) break; // Early exit when target reached

        for (size_t i = 0; i < m_adj[node].size(); ++i) {
            int neighbor = m_adj[node][i];
            double edgeDist = m_adj_distances[node][i];

            if (dist[node] + edgeDist < dist[neighbor]) {
                dist[neighbor] = dist[node] + edgeDist;
                prev[neighbor] = node;
                pq.emplace(dist[neighbor], neighbor);
            }
        }
    }

    // Build the path
    std::vector<int> path;
    path.reserve(n);
    for (int at = v; at != -1; at = prev[at]) {
        path.push_back(at);
    }
    std::reverse(path.begin(), path.end());

    // Optional: check if reachable
    if (path.front() != u) return {}; // empty path if no connection

    return path; // NRVO ensures no copy
}

void Graph::InitializeMemberByParser(int maxNodeIdSeen) {
    IGraph::n = maxNodeIdSeen + 1;
    IGraph::m = 0;
    m_adj.clear();
    m_adj.resize(n);
    m_adj_distances.resize(n);

    m_adj_capacities.resize(n);
    IGraph::vertices.resize(n);
    std::iota(IGraph::vertices.begin(), IGraph::vertices.end(), 0); // fill vertices with 0, 1, ..., maxNodeIdSeen
}

void Graph::print() const {
    for (size_t i = 0; i < m_adj.size(); ++i) {
        std::cout << "Node " << i << ": ";
        for (size_t j = 0; j < m_adj[i].size(); ++j) {
            std::cout << "(" << m_adj[i][j] << ", " << m_adj_capacities[i][j] << ", " << m_adj_distances[i][j] << ") ";
        }
        std::cout << std::endl;
    }
}


bool Graph::IsDistanceMatrixComputed() const {
    return !m_distanceMatrix.empty();
}

void Graph::createDistanceMatrix() {
    int n = getNumNodes();
    const double INF = std::numeric_limits<double>::infinity();

    m_distanceMatrix.assign(n, std::vector<double>(n, INF));
    m_next.assign(n, std::vector<int>(n, -1));

    // --- Initialization ---
    for (int i = 0; i < n; ++i) {
        m_distanceMatrix[i][i] = 0.0;
        m_next[i][i] = i;

        for (size_t k = 0; k < m_adj[i].size(); ++k) {
            int j = m_adj[i][k];
            double dist = m_adj_distances[i][k];
            if (dist < m_distanceMatrix[i][j]) {
                m_distanceMatrix[i][j] = dist;
                m_distanceMatrix[j][i] = dist;   // <== ensure symmetry
                m_next[i][j] = j;
                m_next[j][i] = i;                // <== reverse direction too
            }
        }
    }

    // --- Floyd–Warshall main loop ---
    for (int k = 0; k < n; ++k) {
        for (int i = 0; i < n; ++i) {
            if (m_distanceMatrix[i][k] == INF) continue; // unreachable
            for (int j = 0; j < n; ++j) {
                if (m_distanceMatrix[k][j] == INF) continue;
                double alt = m_distanceMatrix[i][k] + m_distanceMatrix[k][j];
                if (alt < m_distanceMatrix[i][j]) {
                    m_distanceMatrix[i][j] = alt;
                    m_next[i][j] = m_next[i][k]; // first step i→k

                    // also maintain the reverse for undirected graphs
                    m_distanceMatrix[j][i] = alt;
                    m_next[j][i] = m_next[j][k];
                }
            }
        }
    }
}


std::vector<int> Graph::getPrecomputedShortestPath(int u, int v) const {
    std::vector<int> path;
    if (m_next.empty() || m_next[u][v] == -1) return path; // no path
    path.push_back(u);
    while (u != v) {
        u = m_next[u][v];
        if (u == -1) { path.clear(); return path; } // disconnected
        path.push_back(u);
    }
    return path;
}

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
    this->vertices.resize(n);

    // Fill vertices with 0, 1, ..., maxNodeIdSeen
    std::iota(this->vertices.begin(), this->vertices.end(), 0);

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
