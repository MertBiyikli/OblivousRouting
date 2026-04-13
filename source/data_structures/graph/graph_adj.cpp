//
// Created by Mert Biyikli on 20.03.26.
//

#include "../../include/data_structures/graph/graph_adj.h"
#include "../../../include/data_structures/priority_queue.h"

#include <iostream>

IGraph::NeighborRange GraphADJList::neighbors(int node) const {
    return NeighborRange{adjList[node].data(), adjList[node].data() + adjList[node].size()};
}


void GraphADJList::addEdge(int u, int v, double cap, double dist) {
    assert(u >= 0 && u < n);
    assert(v >= 0 && v < n);

    adjList[u].push_back(v);
    capacity[u].push_back(cap);
    distance[u].push_back(dist);

    adjList[v].push_back(u);
    capacity[v].push_back(cap);
    distance[v].push_back(dist);
}


void GraphADJList::finalize() {
    // sort the adjacency lists and capacities and distances based on the neighbor node ids for efficient access
    for (int u = 0; u < n; ++u) {
        std::vector<std::tuple<int, double, double>> neighbors;
        for (size_t i = 0; i < adjList[u].size(); ++i) {
            neighbors.emplace_back(adjList[u][i], capacity[u][i], distance[u][i]);
        }

        std::sort(neighbors.begin(), neighbors.end(),
                  [](const std::tuple<int, double, double>& a, const std::tuple<int, double, double>& b) {
                      return std::get<0>(a) < std::get<0>(b);
                  });

        for (size_t i = 0; i < neighbors.size(); ++i) {
            adjList[u][i] = std::get<0>(neighbors[i]);
            capacity[u][i] = std::get<1>(neighbors[i]);
            distance[u][i] = std::get<2>(neighbors[i]);
        }
    }

    edge_ids.resize(n);
    int edge_id = 0;
    for (int u = 0; u < n; ++u) {
        edge_ids[u].resize(adjList[u].size());
        for (size_t i = 0; i < adjList[u].size(); ++i) {
            edge_ids[u][i] = edge_id++;
        }
    }

    // Build edge_id_start: edge_id_start[u] is the first edge id for node u; last element = total edges
    edge_id_start.resize(n + 1);
    int cur = 0;
    for (int u = 0; u < n; ++u) {
        edge_id_start[u] = cur;
        cur += static_cast<int>(edge_ids[u].size());
    }
    edge_id_start[n] = cur;

    m = edge_id; // total number of directed edges
    is_processed = true;
}

int GraphADJList::getNumDirectedEdges() const {
    return m;
}


const int GraphADJList::getEdgeId(int u, int v) const {
    int edge_id = INVALID_EDGE_ID;

    // Binary search in the sorted adjacency list of u to find v
    const auto& neighbors_u = adjList[u];
    size_t left = 0, right = neighbors_u.size();
    while (left < right) {
        size_t mid = left + (right - left) / 2;
        if (neighbors_u[mid] == v) {
            edge_id = edge_ids[u][mid];
            break;
        } else if (neighbors_u[mid] < v) {
            left = mid + 1;
        } else {
            right = mid;
        }
    }

    return edge_id;
}

const int GraphADJList::getAntiEdge(int e) const {
    const auto& [u, v] = getEdgeEndpoints(e);
    int anti_e = getEdgeId(v, u);
    if (anti_e == INVALID_EDGE_ID) {
        throw std::runtime_error("getAntiEdge: reverse edge not found");
    }
    return anti_e;
}



std::pair<int,int> GraphADJList::getEdgeEndpoints(int e) const {
    if (!is_processed) {
        throw std::runtime_error("edgeEndpoints: graph not finalized");
    }
    if (e < 0 || e >= m) {
        throw std::out_of_range("edgeEndpoints: edge id out of range");
    }

    // Find owning node u via binary search over edge_id_start
    auto it = std::upper_bound(edge_id_start.begin(), edge_id_start.end(), e);
    // upper_bound returns iterator to first element > e; we want previous index
    int idx = static_cast<int>(std::distance(edge_id_start.begin(), it)) - 1;
    if (idx < 0 || idx >= n) {
        throw std::runtime_error("edgeEndpoints: failed to locate owning node");
    }

    int offset = e - edge_id_start[idx];
    // offset should be within edge_ids[idx].size()
    if (offset < 0 || static_cast<size_t>(offset) >= edge_ids[idx].size()) {
        throw std::runtime_error("edgeEndpoints: inconsistent edge id mapping");
    }

    int u = idx;
    int v = adjList[u][offset];
    return {u, v};
}


double GraphADJList::getEdgeCapacity(int u, int v) const {
    double capacity = -1.0;

    // Binary search in the sorted adjacency list of u to find v
    const auto& neighbors_u = adjList[u];
    size_t left = 0, right = neighbors_u.size();
    while (left < right) {
        size_t mid = left + (right - left) / 2;
        if (neighbors_u[mid] == v) {
            capacity = this->capacity[u][mid];
            break;
        } else if (neighbors_u[mid] < v) {
            left = mid + 1;
        } else {
            right = mid;
        }
    }

    if (capacity < 0) {
        throw std::runtime_error("getEdgeCapacity: edge not found");
    }
    return capacity;
}


double GraphADJList::getEdgeDistance(int u, int v) const {
    double distance = -1.0;

    // Binary search in the sorted adjacency list of u to find v
    const auto& neighbors_u = adjList[u];
    size_t left = 0, right = neighbors_u.size();
    while (left < right) {
        size_t mid = left + (right - left) / 2;
        if (neighbors_u[mid] == v) {
            distance = this->distance[u][mid];
            break;
        } else if (neighbors_u[mid] < v) {
            left = mid + 1;
        } else {
            right = mid;
        }
    }

    if (distance < 0) {
        throw std::runtime_error("getEdgeDistance: edge not found");
    }
    return distance;
}

bool GraphADJList::updateEdgeDistance(int u, int v, double dist) {
    bool ok = false;

    const auto& neighbors_u = adjList[u];
    size_t left = 0, right = neighbors_u.size();
    while (left < right) {
        size_t mid = left + (right - left) / 2;
        if (neighbors_u[mid] == v) {
            this->distance[u][mid] = dist;
            ok = true;
            break;
        } else if (neighbors_u[mid] < v) {
            left = mid + 1;
        } else {
            right = mid;
        }
    }

    return ok;
}

bool GraphADJList::updateEdgeCapacity(int u, int v, double cap) {
    bool ok;

    const auto& neighbors_u = adjList[u];
    size_t left = 0, right = neighbors_u.size();
    while (left < right) {
        size_t mid = left + (right - left) / 2;
        if (neighbors_u[mid] == v) {
            this->capacity[u][mid] = cap;
            ok = true;
            break;
        } else if (neighbors_u[mid] < v) {
            left = mid + 1;
        } else {
            right = mid;
        }
    }

    return ok;
}

double GraphADJList::getEdgeCapacity(int e) const {
    auto [u, v] = getEdgeEndpoints(e);
    return getEdgeCapacity(u, v);
}

double GraphADJList::getEdgeDistance(int e) const {
    auto [u, v] = getEdgeEndpoints(e);
    return getEdgeDistance(u, v);
}

bool GraphADJList::updateEdgeDistance(int e, double dist) {
    auto [u, v] = getEdgeEndpoints(e);
    return updateEdgeDistance(u, v, dist);
}


double GraphADJList::getShortestDistance(int u, int v) const {
    double result = std::numeric_limits<double>::infinity();

    // run Dijkstra's algorithm from u to find distance to v
    std::vector<double> dist(n, std::numeric_limits<double>::infinity());
    dist[u] = 0.0;
    MinHeap<double, int> pq(n);
    pq.insert(u, 0.0);

    while (!pq.empty()) {
        int node = pq.top();
        double du = pq.topKey();
        pq.deleteTop();

        if (du > dist[node]) continue; // stale entry
        if (node == v) {
            result = du;
            break;
        }

        for (size_t i = 0; i < adjList[node].size(); ++i) {
            int neighbor = adjList[node][i];
            double weight = distance[node][i];
            double new_dist = du + weight;
            if (new_dist < dist[neighbor]) {
                dist[neighbor] = new_dist;
                pq.insertOrAdjustKey(neighbor, new_dist);
            }
        }
    }
    return result;
}

// maybe think of a way better/faster way of computin the diameter
const double GraphADJList::getDiameter() const {
    double diameter = 0.0;

    // compute the diameter by running Dijkstra's from each node and taking the max distance
    for (int u = 0; u < n; ++u) {
        for (int v = u + 1; v < n; ++v) {
            double dist_uv = getShortestDistance(u, v);
            if (dist_uv > diameter) {
                diameter = dist_uv;
            }
        }
    }

    return diameter;
}

const double GraphADJList::getDiameterApprox() const {
    if (n == 0) return 0.0;
    if (n == 1) return 0.0;

    const double INF = std::numeric_limits<double>::infinity();

    // --- PHASE 1: Two-sweep heuristic for fast approximation ---
    // Start from an arbitrary node (0) and find the farthest node
    std::vector<double> dist(n, INF);
    std::vector<int> parent(n, -1);
    dist[0] = 0.0;

    MinHeap<double, int> pq(n);
    pq.insert(0, 0.0);

    while (!pq.empty()) {
        int node = pq.top();
        double du = pq.topKey();
        pq.deleteTop();

        if (du > dist[node]) continue; // stale entry

        for (size_t i = 0; i < adjList[node].size(); ++i) {
            int neighbor = adjList[node][i];
            double weight = distance[node][i];
            double new_dist = du + weight;
            if (new_dist < dist[neighbor]) {
                dist[neighbor] = new_dist;
                pq.insertOrAdjustKey(neighbor, new_dist);
            }
        }
    }

    // Find farthest node from node 0
    int farthest_u = 0;
    double max_dist_from_0 = 0.0;
    for (int i = 0; i < n; ++i) {
        if (dist[i] < INF && dist[i] > max_dist_from_0) {
            max_dist_from_0 = dist[i];
            farthest_u = i;
        }
    }

    // --- PHASE 2: From the farthest node, find the next farthest node ---
    std::fill(dist.begin(), dist.end(), INF);
    std::fill(parent.begin(), parent.end(), -1);
    dist[farthest_u] = 0.0;

    pq.insert(farthest_u, 0.0);

    while (!pq.empty()) {
        int node = pq.top();
        double du = pq.topKey();
        pq.deleteTop();

        if (du > dist[node]) continue; // stale entry

        for (size_t i = 0; i < adjList[node].size(); ++i) {
            int neighbor = adjList[node][i];
            double weight = distance[node][i];
            double new_dist = du + weight;
            if (new_dist < dist[neighbor]) {
                dist[neighbor] = new_dist;
                pq.insertOrAdjustKey(neighbor, new_dist);
            }
        }
    }

    // Find the diameter as the maximum distance from farthest_u
    double diameter = 0.0;
    for (int i = 0; i < n; ++i) {
        if (dist[i] < INF && dist[i] > diameter) {
            diameter = dist[i];
        }
    }

    return diameter;
}



std::vector<int> GraphADJList::getShortestPath(int u, int v) const {
    // Build the path
    std::vector<int> path;

    // run Dijkstra's algorithm from u to find path to v
    std::vector<double> dist(n, std::numeric_limits<double>::infinity());
    std::vector<int> parent(n, -1);
    dist[u] = 0.0;

    MinHeap<double, int> pq(n);
    pq.insert(u, 0.0);

    while (!pq.empty()) {
        int node = pq.top();
        double du = pq.topKey();
        pq.deleteTop();

        if (du > dist[node]) continue; // stale entry
        if (node == v) {
            break;
        }

        for (size_t i = 0; i < adjList[node].size(); ++i) {
            int neighbor = adjList[node][i];
            double weight = distance[node][i];
            double new_dist = du + weight;
            if (new_dist < dist[neighbor]) {
                dist[neighbor] = new_dist;
                parent[neighbor] = node;
                pq.insertOrAdjustKey(neighbor, new_dist);
            }
        }
    }

    // If v is unreachable, return empty path
    if (dist[v] == std::numeric_limits<double>::infinity()) {
        return path;
    }else {
        // Reconstruct path from v to u using parent pointers
        for (int cur = v; cur != -1; cur = parent[cur]) {
            path.push_back(cur);
        }
        std::reverse(path.begin(), path.end());
    }

    return path; // NRVO ensures no copy
}

void GraphADJList::initializeMemberByParser(int maxNodeIdSeen) {
    IGraph::n = maxNodeIdSeen + 1;
    IGraph::m = 0;

    adjList.resize(n);
    capacity.resize(n);

    distance.resize(n);
    IGraph::vertices.resize(n);
    std::iota(IGraph::vertices.begin(), IGraph::vertices.end(), 0);
}

std::vector<int>
        GraphADJList::getShortestPath(int s, int t, const std::vector<double>& dist_e) const {

    std::vector<int> path;
    // run Dijkstra's algorithm from s to find path to t
    std::vector<double> dist(n, std::numeric_limits<double>::infinity());
    std::vector<int> parent(n, -1);
    dist[s] = 0.0;

    MinHeap<double, int> pq(n);
    pq.insert(s, 0.0);

    while (!pq.empty()) {
        int node = pq.top();
        double du = pq.topKey();
        pq.deleteTop();

        if (du > dist[node]) continue; // stale entry
        if (node == t) {
            break;
        }

        for (size_t i = 0; i < adjList[node].size(); ++i) {
            int neighbor = adjList[node][i];
            double weight = dist_e[edge_ids[node][i]];
            double new_dist = du + weight;
            if (new_dist < dist[neighbor]) {
                dist[neighbor] = new_dist;
                parent[neighbor] = node;
                pq.insertOrAdjustKey(neighbor, new_dist);
            }
        }
    }

    // If t is unreachable, return empty path
    if (dist[t] == std::numeric_limits<double>::infinity()) {
        return path;
    } else {
        // Reconstruct path from t to s using parent pointers
        for (int cur = t; cur != -1; cur = parent[cur]) {
            path.push_back(cur);
        }
        std::reverse(path.begin(), path.end());
    }
    return path;
}


void GraphADJList::print() const {
    for (size_t i = 0; i < adjList.size(); ++i) {
        std::cout << "Node " << i << ": ";
        for (size_t j = 0; j < adjList[i].size(); ++j) {
            std::cout << "(" << adjList[i][j] << ", " << capacity[i][j] << ") ";
        }
        std::cout << std::endl;
    }
}

void GraphADJList::printGraphType() const {
    std::cout << "GraphADJList (Adjacency List) format" << std::endl;
}