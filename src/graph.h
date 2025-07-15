//
// Created by Mert Biyikli on 11.05.25.
//

#ifndef OBLIVOUSROUTING_GRAPH_H
#define OBLIVOUSROUTING_GRAPH_H

#include <vector>
#include <string>
#include <map>
#include <memory>
#include <numeric>
#include <iostream>

class RaeckeGraph{
    int m = 0; // Number of edges
    std::vector<int> m_vertices; // List of vertices
    std::vector<std::vector<int>> m_adj;
    std::vector<std::vector<double>> m_adj_capacities, m_adj_distances;
    std::vector<std::vector<double>> m_distanceMatrix;
public:
    RaeckeGraph() = default;
    RaeckeGraph(int n) : m_adj(n), m_adj_capacities(n), m_adj_distances(n) {
        if(m_vertices.empty()) {
            m_vertices.resize(m_adj.size());
            std::iota(m_vertices.begin(), m_vertices.end(), 0); // Fill vertices with 0, 1, ..., maxNodeIdSeen
        }
    };

    RaeckeGraph(const RaeckeGraph& other)
        : m_vertices(other.m_vertices),
          m_adj(other.m_adj),
          m_adj_capacities(other.m_adj_capacities),
          m_adj_distances(other.m_adj_distances),
          m_distanceMatrix(other.m_distanceMatrix),
          m(other.getNumEdges()){};

    RaeckeGraph& operator=(const RaeckeGraph& other) {
        if (this != &other) {
            m_vertices = other.m_vertices;
            m_adj = other.m_adj;
            m_adj_capacities = other.m_adj_capacities;
            m_adj_distances = other.m_adj_distances;
            m_distanceMatrix = other.m_distanceMatrix;
            m = other.getNumEdges();
        }
        return *this;
    }

    int getNumNodes() const {
        return m_adj.size();
    }

    int getNumEdges() const {
        return m;
    }


    std::vector<int> getVertices() const {
        return m_vertices;
    }

    std::vector<int>& getVertices() {
        return (m_vertices);
    }

    std::vector<double> getDistances(int u) const {
        if(u < 0 || u >= m_adj.size()) {
            throw std::out_of_range("Node index out of range");
        }
        return m_distanceMatrix[u];
    }

    const std::vector<int>& neighbors(int u) const {
        if(u < 0 || u >= m_adj.size()) {
            throw std::out_of_range("Node index out of range");
        }
        return m_adj[u];
    }


    void addAdd(int u , int v, double capacity, double distance = 1.0) {
        if(u >= m_adj.size() || v >= m_adj.size()) {
            throw std::out_of_range("Node index out of range");
        }
        m_adj[u].push_back(v);
        m_adj_capacities[u].push_back(capacity);
        m_adj_distances[u].push_back(distance);

        m_adj[v].push_back(u);
        m_adj_capacities[v].push_back(capacity);
        m_adj_distances[v].push_back(distance);
        m++;
    }

    double getEdgeCapacity(int u, int v) const {
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

    double getEdgeDistance(int u, int v) const {
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

    void updateEdgeDistance(int u, int v, double distance) {
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

    void updateEdgeCapacity(int u, int v, double capacity) {
        if(u >= m_adj.size() || v >= m_adj.size()) {
            throw std::out_of_range("Node index out of range");
        }
        auto it = std::find(m_adj[u].begin(), m_adj[u].end(), v);
        if(it != m_adj[u].end()) {
            int index = std::distance(m_adj[u].begin(), it);
            m_adj_capacities[u][index] = capacity;
            m_adj_capacities[v][index] = capacity; // Assuming undirected graph, update both directions
        } else {
            throw std::runtime_error("Edge not found");
        }
    }

    double getShortestDistance(int u, int v) {
        return m_distanceMatrix[u][v];
    }

    void print() const {
        for (size_t i = 0; i < m_adj.size(); ++i) {
            std::cout << "Node " << i << ": ";
            for (size_t j = 0; j < m_adj[i].size(); ++j) {
                std::cout << "(" << m_adj[i][j] << ", " << m_adj_capacities[i][j] << ", " << m_adj_distances[i][j] << ") ";
            }
            std::cout << std::endl;
        }
    }

    bool IsDistanceMatrixComputed() const {
        return !m_distanceMatrix.empty();
    }

    void createDistanceMatrix() {
        int n = getNumNodes();
        m_distanceMatrix.resize(n, std::vector<double>(n, std::numeric_limits<double>::infinity()));

        for(int i = 0; i < n; ++i) {
            m_distanceMatrix[i][i] = 0.0; // Distance to self is 0
            for(size_t j = 0; j < m_adj[i].size(); ++j) {
                int neighbor = m_adj[i][j];
                double distance = m_adj_distances[i][j];
                m_distanceMatrix[i][neighbor] = distance;
            }
        }

        // Floyd-Warshall algorithm to compute all-pairs shortest paths
        for(int k = 0; k < n; ++k) {
            for(int i = 0; i < n; ++i) {
                for(int j = 0; j < n; ++j) {
                    if(m_distanceMatrix[i][j] > m_distanceMatrix[i][k] + m_distanceMatrix[k][j]) {
                        m_distanceMatrix[i][j] = m_distanceMatrix[i][k] + m_distanceMatrix[k][j];
                    }
                }
            }
        }
    }

    double GetDiameter() const {
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

    std::vector<int> getShortestPath(int u, int v) {
        if (u < 0 || u >= m_adj.size() || v < 0 || v >= m_adj.size()) {
            throw std::out_of_range("Node index out of range");
        }

        std::vector<int> path;

        // compute shortest path using dijskta
        std::vector<double> dist(m_adj.size(), std::numeric_limits<double>::infinity());
        std::vector<int> prev(m_adj.size(), -1);
        dist[u] = 0.0;
        std::vector<bool> visited(m_adj.size(), false);
        for (int i = 0; i < m_adj.size(); ++i) {
            int minNode = -1;
            double minDist = std::numeric_limits<double>::infinity();
            for (int j = 0; j < m_adj.size(); ++j) {
                if (!visited[j] && dist[j] < minDist) {
                    minDist = dist[j];
                    minNode = j;
                }
            }

            if (minNode == -1) break; // No more reachable nodes

            visited[minNode] = true;

            for (size_t j = 0; j < m_adj[minNode].size(); ++j) {
                int neighbor = m_adj[minNode][j];
                double edgeDistance = m_adj_distances[minNode][j];
                if (dist[minNode] + edgeDistance < dist[neighbor]) {
                    dist[neighbor] = dist[minNode] + edgeDistance;
                    prev[neighbor] = minNode;
                }
            }
        }

        // return the path from u to v
        for (int at = v; at != -1; at = prev[at]) {
            path.push_back(at);
        }
        return path;
    }

};

struct Arc{
    int source;
    int target;
    int id;
    int rev_id;
    double capacity;

    int GetId() const {return id;};
    int GetReverseArcId() const {return rev_id;};
    int GetSource() const {return source;};
    int GetTarget() const {return target;};
    double GetCapacity() const {return capacity;};
    bool operator!=(const Arc& arc){return (this->GetId() == arc.GetId());};
};

struct Edge{
    int source;
    int target;
    double capacity;
    double distance = 0.0; // optional distance, used in some algorithms
};

class DiGraph;

class Graph{
public:
    Graph() : m_iNumNodes(0), m_iNumEdges(0) {};
    Graph(int n): m_iNumNodes(n), m_iNumEdges(0), m_adj(n), m_vertices(n) {std::iota(std::begin(m_vertices), std::end(m_vertices), 0);}
    Graph(const Graph& graph) {
        m_iNumNodes = graph.m_iNumNodes;
        m_iNumEdges = graph.m_iNumEdges;
        m_adj = graph.m_adj;
        m_distanceMatrix = graph.m_distanceMatrix;
        m_ArcMatrix = graph.m_ArcMatrix;
        m_vertices = graph.m_vertices;
        if(graph.m_diGraph) {
            m_diGraph = std::make_unique<DiGraph>(*graph.m_diGraph);
        } else {
            m_diGraph.reset();
        }
    }

    RaeckeGraph getRaeckeGraph() const {
        RaeckeGraph rg(this->numNodes());
        for(int i = 0; i < m_iNumNodes; ++i) {
            for(const auto& edge : m_adj[i]) {
                if(edge->source >= edge->target) continue; // Avoid adding the same undirected edge twice
                rg.addAdd(edge->source, edge->target, edge->capacity);
            }
        }
        return rg;
    }


    void addEdge(int source, int target, double capacity);


    const std::vector<std::shared_ptr<Edge>>& neighbors(int node) const {
        return m_adj[node];
    }

    std::vector<std::shared_ptr<Edge>>& neighbors(int node) {
        return (m_adj[node]);
    }

    int numNodes() const;
    int numEdges() const;

    const std::vector<int>& GetVertices() const;
    std::vector<int>& GetVertices();


    std::shared_ptr<Edge> GetReverseEdge(const std::shared_ptr<Edge>& edge) {
        if(edge) {
            m_adj[edge->target];
            for(auto& e : m_adj[edge->target]) {
                if(e->target == edge->source && e->source == edge->target) {
                    return e;
                }
            }
        }
        throw std::runtime_error("Reverse edge with id " + std::to_string(edge->target) + " " + std::to_string(edge->source) + " not found.");
    }

    // Todo
    // bool isValid() const;
    void readGraph(const std::string& filename);

    void readLFGFile(const std::string& filename, bool withDistances);
    void print() const;
    const DiGraph& GetDiGraph() const; // TODO: implement directed graph


    void createDistanceMatrix();
    const std::vector<std::vector< double> >& GetDistanceMatrix() const;
    const std::unordered_map<int, std::shared_ptr<Edge>>& GetArcMatrix(int node) const;
    double GetGraphDiameter() const;

    std::vector<std::shared_ptr<Edge>> getShortestPath(int source, int target) const;

    std::vector<std::vector< double> > m_distanceMatrix;
    std::vector<std::unordered_map< int, std::shared_ptr<Edge>> > m_ArcMatrix;



private:
    std::vector<std::vector<std::shared_ptr<Edge> > > m_adj;
    std::vector<int> m_vertices;

    // directed symmetric graph representation
    mutable std::unique_ptr<DiGraph> m_diGraph;

    int m_iNumNodes;
    int m_iNumEdges;

};

class DiGraph{
public:
    DiGraph() : m_iNumNodes(0), m_iNumArcs(0) {};
    DiGraph(int n): m_iNumNodes(n), m_iNumArcs(0)/** TODO ***/ {}


    void addUndirectedArc(int source, int target, double capacity);


    const std::vector<Arc>& outgoingNeighbors(int node) const;

    int numNodes() const;
    int numArcs() const;

    bool isValid();
    bool HasArcId(int id) const;
    const Arc& GetArcById(int id) const;

    void readDiGraph(const std::string& filename);
    void print() const;

    const std::unordered_map<int, Arc>& GetArcs() const;

private:
    //std::vector<std::vector<Edge> > m_adj;

    // important as map, since we need to represent undirected edges as directed edges(arcs)
    std::unordered_map<int, Arc> m_IdToArc;

    int m_idArc = 0;
    int m_iNumNodes;
    int m_iNumArcs;
};


#endif //OBLIVOUSROUTING_GRAPH_H
