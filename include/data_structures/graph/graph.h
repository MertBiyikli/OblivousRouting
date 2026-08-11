//
// Created by Mert Biyikli on 27.07.26.
//

#ifndef OBLIVIOUSROUTING_GRAPH_H
#define OBLIVIOUSROUTING_GRAPH_H
#include <algorithm>
#include <cassert>
#include <iostream>
#include <numeric>
#include <queue>
#include <span>
#include <stack>
#include <vector>

#include "core/errors.h"
#include "../../../include/data_structures/priority_queue.h"


struct EdgeData {
    double capacity;
    double weight;
    friend std::ostream& operator<<(std::ostream& os, const EdgeData& obj) {
        os << "cap: "<<  obj.capacity << " weigh: " << obj.weight << "\n";
        return os;
    }

    void resetEdgeWeight(double _new_weight = 1.0) {
        weight = _new_weight;
    }
};

namespace  optimized {
    // Lightweight entry stored in an adjacency list.
    struct Edge {
        int tail;
        int id;
    };


    enum class UndoKind{
        edge,
        node,
        data
    };

    struct Undo {
        UndoKind kind;
        int padding[3]{};
        int id;
    };


    template<typename Data>
    class Graph {
    public:
        struct ReindexedSubgraph;

        struct InputEdge {
            int head;
            int tail;
            Data data;
        };
        struct SubgraphLevel {
            // State before selecting the induced subgraph.
            int parent;

            // Baseline of the selected induced subgraph. Normal rollback is not
            // allowed to cross this point; restoreSubgraph() deliberately does.
            int baseline;
        };
    private:

        struct StoredEdge {
            int head;
            int tail;
            Data data;
            int headIndex;
            int tailIndex;
            bool alive = true;
        };


        std::vector<std::vector<Edge>> edges;
        std::vector<int> edgeBound;

        std::vector<int> nodes;
        std::vector<int> nodeIndex;
        int nodeBound = 0;

        // One record per undirected edge. The ID never changes.
        std::vector<StoredEdge> storedEdges;
        int activeEdgeCount = 0;

        std::vector<Undo> undoLog;
        std::vector<Data> oldData;
        std::vector<SubgraphLevel> subgraphLevels;

        // Reused by Dijkstra and by the path-edge reconstruction helpers.
        mutable std::vector<double> distanceBuffer;
        mutable std::vector<int> parentBuffer;
        mutable std::vector<int> parentEdgeBuffer;
        mutable int parentSource = -1;
        mutable std::uint64_t parentGraphVersion = 0;
        std::uint64_t graphVersion = 0;

    public:
        Graph(int n, const std::vector<InputEdge>& input)
       : edges(n),
         edgeBound(n),
         nodes(n),
         nodeIndex(n),
         nodeBound(n) {
            std::iota(nodes.begin(), nodes.end(), 0);
            std::iota(nodeIndex.begin(), nodeIndex.end(), 0);
            storedEdges.reserve(input.size());

            for (const auto& edge : input) {


                if (edge.head == edge.tail) {
                    throw std::invalid_argument("Self-loops are not supported.");
                }

                const int id = static_cast<int>(storedEdges.size());
                const int headIndex = edges[edge.head].size();
                const int tailIndex = edges[edge.tail].size();

                edges[edge.head].push_back({edge.tail, 2 * id});
                edges[edge.tail].push_back({edge.head, 2 * id + 1});

                storedEdges.push_back({
                    edge.head,
                    edge.tail,
                    edge.data,
                    headIndex,
                    tailIndex
                });
            }

            for (std::size_t u = 0; u < n; ++u) {
                edgeBound[u] = edges[u].size();
            }

            activeEdgeCount = storedEdges.size();
        }

        const int getNumNodes() const{
            return nodeBound;
        }

        int getNumNodes(){
            return nodeBound;
        }

        // Same iteration style as the posted Graph class.
        auto begin() noexcept { return nodes.begin(); }
        auto end() noexcept { return nodes.begin() + nodeBound; }
        auto begin() const noexcept { return nodes.cbegin(); }
        auto end() const noexcept { return nodes.cbegin() + nodeBound; }
        auto cbegin() const noexcept { return begin(); }
        auto cend() const noexcept { return end(); }

        auto beginEdge(int u) noexcept {
            return edges[u].begin();
        }

        auto endEdge(int u) noexcept {
            return edges[u].begin() + edgeBound[u];
        }

        auto cbeginEdge(int u) const noexcept {
            return edges[u].cbegin();
        }

        auto cendEdge(int u) const noexcept {
            return edges[u].cbegin() + edgeBound[u];
        }

        std::span<Edge> edgesOf(int u) noexcept {
            return {beginEdge(u), endEdge(u)};
        }

        std::span<const Edge> edgesOf(int u) const noexcept {
            return {cbeginEdge(u), cendEdge(u)};
        }

        Edge& getEdge(int u, int index) noexcept {
            assert(index < edgeBound[u]);
            return edges[u][index];
        }

        const Edge& getEdge(int u, int index) const noexcept {
            assert(index < edgeBound[u]);
            return edges[u][index];
        }

        Edge& reverse(const Edge& edge) noexcept {
            StoredEdge& stored = storedEdges[toUndirectedEdgeId(edge.id)];
            return edge.tail == stored.tail
                ? edges[stored.tail][stored.tailIndex]
                : edges[stored.head][stored.headIndex];
        }

        const Edge& reverse(const Edge& edge) const noexcept {
            const StoredEdge& stored = storedEdges[toUndirectedEdgeId(edge.id)];
            return edge.tail == stored.tail
                ? edges[stored.tail][stored.tailIndex]
                : edges[stored.head][stored.headIndex];
        }

        const Edge reverse(const int& edge) const noexcept {
            const StoredEdge& stored = storedEdges[toUndirectedEdgeId(edge)];
            if ((edge & 1) == 0) {
                return edges[stored.tail][stored.tailIndex];
            }
            return edges[stored.head][stored.headIndex];
        }

        int size() const noexcept {
            return nodeBound;
        }

        int getNumDirectedEdges() const noexcept {
            return activeEdgeCount*2;
        }

        int getNumUndirectedEdges() const noexcept {
            return static_cast<int>(activeEdgeCount);
        }

        int globalEdgeCount() const noexcept {
            return storedEdges.size();
        }

        bool alive(int u) const noexcept {
            return u >= 0
                && static_cast<int>(u) < nodeIndex.size()
                && nodeIndex[u] < nodeBound;
        }

        bool edgeAlive(int id) const noexcept {
            if (id < 0) {
                return false;
            }
            const int undirected = toUndirectedEdgeId(id);
            return undirected < static_cast<int>(storedEdges.size()) &&
                storedEdges[undirected].alive;
        }

        int degree(int u) const noexcept {
            return edgeBound[u];
        }

        int globalDegree(int u) const noexcept {
            return edges[u].size();
        }

        const Data& edgeData(int id) const {
            return getStoredEdge(id).data;
        }

        Data& edgeData(int id) {
            return getStoredEdge(id).data;
        }

        std::pair<int, int> getEdgeEndpoints(int id) const {
            const auto& edge = getStoredEdge(id);
            if ((id & 1) == 0) {
                return {edge.head, edge.tail};
            }
            return {edge.tail, edge.head};
        }

        // O(1)
        bool removeEdge(int id) {
            if (!edgeAlive(id)) {
                return false;
            }

            const int undirected = toUndirectedEdgeId(id);
            removeEdgeWithoutLogging(undirected);
            undoLog.push_back({UndoKind::edge, {}, undirected});
            return true;
        }

        // O(active degree of u)
        bool remove(int u) {
            if (!alive(u)) {
                return false;
            }

            while (edgeBound[u] > 0) {
                removeEdge(edges[u][edgeBound[u] - 1].id);
            }

            const int head = nodeIndex[u];
            const int tail = nodeBound - 1;

            swapnodes(head, tail);
            --nodeBound;

            // Log after the edges. Rollback therefore restores the nodex first.
            undoLog.push_back({
                UndoKind::node,
                {},
                static_cast<int>(u)
            });

            return true;
        }

        /**
     * Focus this Graph on the subgraph induced by the given active vertices.
     *
     * Vertex and edge IDs remain unchanged. No second Graph is allocated:
     * vertices outside the subset and their incident edges are moved behind
     * the existing active bounds.
     *
     * The subset must contain distinct vertices that are active in the current
     * graph. Calls may be nested and are undone one level at a time with
     * restoreSubgraph().
     *
     * Complexity: O(n + sum of degrees of removed vertices).
     */
    template<std::input_iterator It>
    Result<void> subgraph(It subsetBegin, It subsetEnd) {
        std::vector<std::uint8_t> selected(nodes.size(), 0);

        // Validate the complete input before changing the graph.
        for (auto it = subsetBegin; it != subsetEnd; ++it) {
            const int u = static_cast<int>(*it);

            if (!alive(u)) {
                return makeErrorMessage(ErrorCode::InvalidArgument, "Subgraph contains an inactive vertex.");
            }

            if (selected[u] != 0) {
                return makeErrorMessage(ErrorCode::InvalidArgument, "Subgraph contains a duplicate vertex.");
            }

            selected[u] = 1;
        }

        const int parent = checkpoint();

        // Reserve before mutating so logging cannot reallocate mid-operation.
        undoLog.reserve(
            undoLog.size() + activeEdgeCount + nodeBound
        );
        subgraphLevels.push_back({parent, parent});

        /*
         * remove() swaps another active vertex into the current position.
         * Therefore the index advances only when the current vertex is kept.
         */
        std::size_t index = 0;
        while (index < nodeBound) {
            const int u = nodes[index];

            if (selected[u] != 0) {
                ++index;
            } else {
                remove(u);
            }
        }

        subgraphLevels.back().baseline = checkpoint();

        checkGraphIntegrity();
            return {};
    }

    Result<void> subgraph(std::span<const int> subset) {
        auto sub = subgraph(subset.begin(), subset.end());
            if (!sub) {
                return getError(sub);
            }
            return {};
    }

        /**
    * Build an independent copy of the active induced subgraph.
    *
    * The returned graph uses contiguous vertex IDs [0, subset size) in the
    * order supplied by the caller. Its edge IDs are also contiguous and are
    * assigned in the order of the corresponding original edge IDs.
    *
    * Only currently active vertices and edges are considered. The subset
    * must contain distinct vertices that are active in this graph.
    */
template<std::input_iterator It>
Result<Graph<Data>::ReindexedSubgraph> reindexedSubgraph(It subsetBegin,It subsetEnd) const {
    std::vector<int> newToOldVertex;
    std::vector<int> oldToNewVertex(nodes.size(),{-1});

    for (auto it = subsetBegin; it != subsetEnd; ++it) {
        const int oldVertex = static_cast<int>(*it);

        if (!alive(oldVertex)) {
            return makeErrorMessage(ErrorCode::InvalidArgument, "Subgraph contains an inactive vertex.");
        }

        if (oldToNewVertex[oldVertex] != -1) {
            return makeErrorMessage(ErrorCode::InvalidArgument, "Subgraph contains a duplicate vertex.");
        }

        const int newVertex = static_cast<int>(newToOldVertex.size());

        oldToNewVertex[oldVertex] = newVertex;
        newToOldVertex.push_back(oldVertex);
    }

    std::vector<InputEdge> inputEdges;
    std::vector<int> newToOldEdge;

    inputEdges.reserve(activeEdgeCount);
    newToOldEdge.reserve(activeEdgeCount);

    for (std::size_t rawId = 0;
         rawId < storedEdges.size();
         ++rawId) {
        const auto& edge = storedEdges[rawId];

        if (!edge.alive) {
            continue;
        }

        const int newFrom = oldToNewVertex[edge.head];
        const int newTo = oldToNewVertex[edge.tail];

        if (newFrom == -1 || newTo == -1) {
            continue;
        }

        inputEdges.push_back({
            newFrom,
            newTo,
            edge.data
        });

        newToOldEdge.push_back(
            static_cast<int>(rawId)
        );
    }

    return ReindexedSubgraph{
        optimized::Graph<Data>(newToOldVertex.size(), inputEdges),
        std::move(newToOldVertex),
        std::move(oldToNewVertex),
        std::move(newToOldEdge)
    };
}

Result<ReindexedSubgraph> reindexedSubgraph(std::span<const int> subset) const {
    return reindexedSubgraph(subset.begin(),subset.end());
}

    /**
     * Restore the graph state that existed immediately before the most recent
     * subgraph() call.
     */
    Result<void> restoreSubgraph() {
        if (subgraphLevels.empty()) {
            return makeErrorMessage(ErrorCode::LogicError, "No subgraph exists to restore.");
        }

        const int parent = subgraphLevels.back().parent;

        rollbackUnchecked(parent);
        subgraphLevels.pop_back();
        checkGraphIntegrity();
            return {};
    }

    int subgraphDepth() const noexcept {
        return subgraphLevels.size();
    }

        void setEdgeData(int id, Data data) {
            const int undirected = toUndirectedEdgeId(id);
            auto& edge = getStoredEdgeByUndirected(undirected);

            oldData.push_back(edge.data);
            undoLog.push_back({UndoKind::data, {}, undirected});
            edge.data = std::move(data);
        }

        int checkpoint() const noexcept {
            return undoLog.size();
        }

        Result<void> rollback(int checkpoint) {
            if (checkpoint > undoLog.size()) {
                return makeErrorMessage(ErrorCode::InvalidArgument, "Invalid graph checkpoint.");
            }

            if (!subgraphLevels.empty() &&
                checkpoint < subgraphLevels.back().baseline) {
                return makeErrorMessage(ErrorCode::InvalidArgument, "Checkpoint belongs to a parent subgraph.");
            }

            rollbackUnchecked(checkpoint);
            return {};
        }

        Result<void> restoreAll() {
            auto roll = rollback(
                subgraphLevels.empty()
                    ? int{0}
                    : subgraphLevels.back().baseline
            );
            if (!roll) {
                return getError(roll);
            }else {
                return {};
            }
        }

        // Forget rollback history and make the current state the new baseline.
        Result<void> commit() {
            if (!subgraphLevels.empty()) {
                return makeErrorMessage(ErrorCode::LogicError, "Cannot commit while a subgraph is active.");
            }

            undoLog.clear();
            oldData.clear();
            return {};
        }

        void checkGraphIntegrity() const {
#ifndef NDEBUG
            int count = 0;

            for (int id = 0; id < storedEdges.size(); ++id) {
                const auto& edge = storedEdges[id];
                const bool activeFrom = edge.headIndex < edgeBound[edge.head];
                const bool activeTo = edge.tailIndex < edgeBound[edge.tail];

                assert(activeFrom == edge.alive);
                assert(activeTo == edge.alive);
                assert(edges[edge.head][edge.headIndex].id == 2 * id);
                assert(edges[edge.tail][edge.tailIndex].id == 2 * id + 1);

                if (edge.alive) {
                    assert(alive(edge.head));
                    assert(alive(edge.tail));
                    ++count;
                }
            }

            assert(count == activeEdgeCount);
#endif
        }
        void print() {
            for (const auto node : (*this)) {
                for (const auto& edge : (*this).edgesOf(node)) {
                    const auto neighbour = edge.tail;
                    const auto edgeId = edge.id;
                    const auto& data = (*this).edgeData(edgeId);

                    std::cout
                        << node << " -> " << neighbour
                        << ", capacity = " << data.capacity
                        << '\n';
                }
            }
        }

        int edgeId(int head, int tail) const {
            assert(head <= nodeBound && tail <= nodeBound && "node is outside of the graph.");
            auto it = std::find_if(storedEdges.begin(), storedEdges.end(), [&](const StoredEdge& edge) {return (edge.head == head &&  edge.tail == tail );});
            if (it != storedEdges.end()) {
                const int undirected = static_cast<int>(std::distance(storedEdges.begin(), it));
                return 2 * undirected;
            }
            it = std::find_if(storedEdges.begin(), storedEdges.end(), [&](const StoredEdge& edge) {return (edge.head == tail &&  edge.tail == head );});
            if (it == storedEdges.end()) {
                return -1;
            }
            const int undirected = static_cast<int>(std::distance(storedEdges.begin(), it));
            return 2 * undirected + 1;
        }

        /**
         * Dijkstra using Data::weight. Node and edge IDs remain the original
         * stable IDs, including while an induced subgraph is active.
         */
        Result<std::vector<int>> getShortestPath(int src, int tgt) const {
            auto valid = validatePathQuery(src, tgt, nullptr);
            if (!valid) {
                return getError(valid);
            }

            runDijkstra(src, tgt, nullptr);
            return reconstructNodePath(src, tgt);
        }

        /**
         * Dijkstra using an external weight vector indexed by stable edge ID.
         */
        Result<std::vector<int>> getShortestPath(int src,int tgt,const std::vector<double>& distance) const {
            auto valid = validatePathQuery(src, tgt, &distance);
            if (!valid) {
                return getError(valid);
            }

            runDijkstra(src, tgt, &distance);
            return reconstructNodePath(src, tgt);
        }

        /**
         * Reconstruct edge IDs after the most recent unidirectional Dijkstra.
         */
        Result<std::vector<int>> getPathEdgesFromParent(int src, int tgt) const {
            auto valid = validatePathQuery(src, tgt, nullptr);
            if (!valid) {
                return getError(valid);
            }

            if (parentSource != src ||
                parentGraphVersion != graphVersion ||
                parentBuffer.size() != nodes.size() ||
                parentEdgeBuffer.size() != nodes.size()) {
                return makeErrorMessage(
                    ErrorCode::LogicError,
                    "No shortest-path search from this source is available."
                );
            }

            std::vector<int> pathEdges;
            if (src == tgt) {
                return pathEdges;
            }

            if (parentBuffer[tgt] == -1) {
                return pathEdges;
            }

            for (int v = tgt; v != src; v = parentBuffer[v]) {
                if (v == -1 || parentEdgeBuffer[v] == -1) {
                    pathEdges.clear();
                    return pathEdges;
                }
                pathEdges.push_back(parentEdgeBuffer[v]);
            }

            std::reverse(pathEdges.begin(), pathEdges.end());
            return pathEdges;
        }

        Result<std::vector<int>> getPathEdges(int src, int tgt) const {
            auto path = getShortestPath(src, tgt);
            if (!path) {
                return getError(path);
            }
            return getPathEdgesFromParent(src, tgt);
        }

        Result<std::vector<int>> getPathEdges(int src,int tgt,const std::vector<double>& distance) const {
            auto path = getShortestPath(src, tgt, distance);
            if (!path) {
                return getError(path);
            }
            return getPathEdgesFromParent(src, tgt);
        }

        Result<double> getShortestDistance(int src, int tgt) const {
            auto valid = validatePathQuery(src, tgt, nullptr);
            if (!valid) {
                return getError(valid);
            }

            runDijkstra(src, tgt, nullptr);
            return distanceBuffer[tgt];
        }

        Result<double> getShortestDistance(int src,int tgt,const std::vector<double>& distance) const {
            auto valid = validatePathQuery(src, tgt, &distance);
            if (!valid) {
                return getError(valid);
            }

            runDijkstra(src, tgt, &distance);
            return distanceBuffer[tgt];
        }

        /**
         * Exact weighted diameter: maximum finite shortest-path distance over
         * all active node pairs. For a disconnected graph this is the largest
         * finite diameter among its active components.
         */
        double getDiameter() const {
            double diameter = 0.0;

            for (const int src : *this) {
                runDijkstra(src, -1, nullptr);
                for (const int tgt : *this) {
                    if (distanceBuffer[tgt] <
                        std::numeric_limits<double>::infinity()) {
                        diameter = std::max(diameter, distanceBuffer[tgt]);
                    }
                }
            }

            return diameter;
        }

        /**
         * Two-sweep weighted-diameter approximation. It is exact on trees and
         * is a lower bound for a connected undirected weighted graph.
         */
        double getDiameterApprox() const {
            if (nodeBound < 2) {
                return 0.0;
            }

            const int start = nodes.front();
            runDijkstra(start, -1, nullptr);
            const int first = farthestReachableNode();

            runDijkstra(first, -1, nullptr);
            const int second = farthestReachableNode();
            return distanceBuffer[second];
        }

        Result<std::vector<int>> getShortestPathBidirectionalSearch(int src,int tgt) const {
            return getShortestPathBidirectionalSearchImpl(src, tgt, nullptr);
        }

        Result<std::vector<int>> getShortestPathBidirectionalSearch(int src,int tgt,const std::vector<double>& distance) const {
            return getShortestPathBidirectionalSearchImpl(src,tgt,&distance);
        }

    private:
        using QueueEntry = std::pair<double, int>;
        using MinQueue = std::priority_queue<QueueEntry,std::vector<QueueEntry>,std::greater<QueueEntry>>;

        Result<void> validatePathQuery( int src,int tgt,const std::vector<double>* distance) const {
            if (!alive(src) || !alive(tgt)) {
                return makeErrorMessage(ErrorCode::InvalidArgument,"Shortest-path endpoint is outside the active graph.");
            }

            if (distance != nullptr &&
                distance->size() < storedEdges.size()) {
                return makeErrorMessage(ErrorCode::InvalidArgument,"Distance vector must contain one value per stable edge ID.");
            }

            return {};
        }

        double edgeWeight(int edgeId,const std::vector<double>* distance) const {
            if (distance != nullptr) {
                return (*distance)[edgeId];
            }
            return static_cast<double>(storedEdges[toUndirectedEdgeId(edgeId)].data.weight);
        }

        void runDijkstra(int src,int tgt,const std::vector<double>* distance) const {
            const double infinity = std::numeric_limits<double>::infinity();
            const std::size_t n = nodes.size();

            distanceBuffer.assign(n, infinity);
            parentBuffer.assign(n, -1);
            parentEdgeBuffer.assign(n, -1);
            parentSource = src;
            parentGraphVersion = graphVersion;

            MinQueue queue;
            distanceBuffer[src] = 0.0;
            queue.emplace(0.0, src);

            while (!queue.empty()) {
                const auto [distanceToU, u] = queue.top();
                queue.pop();

                if (distanceToU > distanceBuffer[u]) {
                    continue;
                }
                if (u == tgt) {
                    break;
                }

                for (const Edge& edge : edgesOf(u)) {
                    const double weight = edgeWeight(edge.id, distance);
                    if (weight < 0.0) {
                        continue;
                    }

                    const int v = edge.tail;
                    const double candidate = distanceToU + weight;
                    if (candidate < distanceBuffer[v]) {
                        distanceBuffer[v] = candidate;
                        parentBuffer[v] = u;
                        parentEdgeBuffer[v] = edge.id;
                        queue.emplace(candidate, v);
                    }
                }
            }
        }

        std::vector<int> reconstructNodePath(int src, int tgt) const {
            const double infinity =
                std::numeric_limits<double>::infinity();
            std::vector<int> path;

            if (distanceBuffer[tgt] == infinity) {
                return path;
            }

            for (int v = tgt; v != -1; v = parentBuffer[v]) {
                path.push_back(v);
            }
            std::reverse(path.begin(), path.end());

            if (path.empty() || path.front() != src) {
                path.clear();
            }
            return path;
        }

        int farthestReachableNode() const {
            int farthest = nodes.front();
            double farthestDistance = -1.0;
            const double infinity =
                std::numeric_limits<double>::infinity();

            for (const int node : *this) {
                if (distanceBuffer[node] < infinity &&
                    distanceBuffer[node] > farthestDistance) {
                    farthest = node;
                    farthestDistance = distanceBuffer[node];
                }
            }
            return farthest;
        }

        Result<std::vector<int>> getShortestPathBidirectionalSearchImpl(int src,int tgt,const std::vector<double>* distance) const {
            auto valid = validatePathQuery(src, tgt, distance);
            if (!valid) {
                return getError(valid);
            }

            if (src == tgt) {
                return std::vector<int>{src};
            }

            const double infinity = std::numeric_limits<double>::infinity();
            const std::size_t n = nodes.size();

            std::vector<double> forwardDistance(n, infinity);
            std::vector<double> backwardDistance(n, infinity);
            std::vector<int> forwardParent(n, -1);
            std::vector<int> backwardParent(n, -1);
            std::vector<std::uint8_t> forwardSettled(n, 0);
            std::vector<std::uint8_t> backwardSettled(n, 0);

            MinQueue forwardQueue;
            MinQueue backwardQueue;
            forwardDistance[src] = 0.0;
            backwardDistance[tgt] = 0.0;
            forwardQueue.emplace(0.0, src);
            backwardQueue.emplace(0.0, tgt);

            int meetingNode = -1;
            double bestPathLength = infinity;

            while (!forwardQueue.empty() && !backwardQueue.empty()) {
                const double forwardMinimum = forwardQueue.top().first;
                const double backwardMinimum = backwardQueue.top().first;
                if (forwardMinimum + backwardMinimum >= bestPathLength) {
                    break;
                }

                if (forwardMinimum <= backwardMinimum) {
                    const auto [distanceToU, u] = forwardQueue.top();
                    forwardQueue.pop();

                    if (distanceToU > forwardDistance[u] ||
                        forwardSettled[u] != 0) {
                        continue;
                    }
                    forwardSettled[u] = 1;

                    if (backwardDistance[u] < infinity) {
                        const double candidate =
                            distanceToU + backwardDistance[u];
                        if (candidate < bestPathLength) {
                            bestPathLength = candidate;
                            meetingNode = u;
                        }
                    }

                    for (const Edge& edge : edgesOf(u)) {
                        const double weight =
                            edgeWeight(edge.id, distance);
                        if (weight < 0.0) {
                            continue;
                        }

                        const int v = edge.tail;
                        const double candidate = distanceToU + weight;
                        if (candidate < forwardDistance[v]) {
                            forwardDistance[v] = candidate;
                            forwardParent[v] = u;
                            forwardQueue.emplace(candidate, v);
                        }

                        if (backwardDistance[v] < infinity) {
                            const double complete =
                                forwardDistance[v] + backwardDistance[v];
                            if (complete < bestPathLength) {
                                bestPathLength = complete;
                                meetingNode = v;
                            }
                        }
                    }
                } else {
                    const auto [distanceToU, u] = backwardQueue.top();
                    backwardQueue.pop();

                    if (distanceToU > backwardDistance[u] ||
                        backwardSettled[u] != 0) {
                        continue;
                    }
                    backwardSettled[u] = 1;

                    if (forwardDistance[u] < infinity) {
                        const double candidate =
                            forwardDistance[u] + distanceToU;
                        if (candidate < bestPathLength) {
                            bestPathLength = candidate;
                            meetingNode = u;
                        }
                    }

                    for (const Edge& edge : edgesOf(u)) {
                        const double weight =
                            edgeWeight(edge.id, distance);
                        if (weight < 0.0) {
                            continue;
                        }

                        const int v = edge.tail;
                        const double candidate = distanceToU + weight;
                        if (candidate < backwardDistance[v]) {
                            backwardDistance[v] = candidate;
                            backwardParent[v] = u;
                            backwardQueue.emplace(candidate, v);
                        }

                        if (forwardDistance[v] < infinity) {
                            const double complete =
                                forwardDistance[v] + backwardDistance[v];
                            if (complete < bestPathLength) {
                                bestPathLength = complete;
                                meetingNode = v;
                            }
                        }
                    }
                }
            }

            std::vector<int> path;
            if (meetingNode == -1) {
                return path;
            }

            for (int v = meetingNode; v != -1; v = forwardParent[v]) {
                path.push_back(v);
            }
            std::reverse(path.begin(), path.end());

            for (int v = backwardParent[meetingNode];
                 v != -1;
                 v = backwardParent[v]) {
                path.push_back(v);
            }

            if (path.empty() || path.front() != src ||
                path.back() != tgt) {
                path.clear();
            }
            return path;
        }


    private:
        Result<void> checknode(int u) const {
            if (u < 0 || static_cast<int>(u) >= nodes.size()) {
                return makeErrorMessage(ErrorCode::InvalidArgument, "node is outside the graph.");
            }
            return {};
        }

        static int toUndirectedEdgeId(int directedId) noexcept {
            return directedId >> 1;
        }

        StoredEdge& getStoredEdgeByUndirected(int id) {
            assert(id >= 0 && id < static_cast<int>(storedEdges.size()) && "Edge is outside the graph.");
            return storedEdges[id];
        }

        const StoredEdge& getStoredEdgeByUndirected(int id) const {
            assert(id >= 0 && id < static_cast<int>(storedEdges.size()) && "Edge is outside the graph.");
            return storedEdges[id];
        }

        StoredEdge& getStoredEdge(int id) {
            return getStoredEdgeByUndirected(toUndirectedEdgeId(id));
        }

        const StoredEdge& getStoredEdge(int id) const {
            return getStoredEdgeByUndirected(toUndirectedEdgeId(id));
        }

        void updateIndex(int owner, int index) noexcept {
            const Edge& adjacencyEdge = edges[owner][index];
            StoredEdge& storedEdge = storedEdges[toUndirectedEdgeId(adjacencyEdge.id)];

            if (storedEdge.head == owner) {
                storedEdge.headIndex = index;
            } else {
                storedEdge.tailIndex = index;
            }
        }

        void swapAdjacencyEdges(int owner,int first,int second) noexcept {
            if (first == second) {
                return;
            }

            std::swap(edges[owner][first], edges[owner][second]);
            updateIndex(owner, first);
            updateIndex(owner, second);
        }

        void removeEdgeWithoutLogging(int id) noexcept {
            StoredEdge& edge = storedEdges[id];

            const int lastFrom = edgeBound[edge.head] - 1;
            swapAdjacencyEdges(edge.head, edge.headIndex, lastFrom);
            --edgeBound[edge.head];

            // edge.tailIndex may have changed during the first swap.
            const int lastTo = edgeBound[edge.tail] - 1;
            swapAdjacencyEdges(edge.tail, edge.tailIndex, lastTo);
            --edgeBound[edge.tail];

            edge.alive = false;
            --activeEdgeCount;
        }

        void restoreEdge(int id) noexcept {
            StoredEdge& edge = storedEdges[id];
            assert(!edge.alive);
            assert(alive(edge.head) && alive(edge.tail));

            swapAdjacencyEdges(edge.head, edge.headIndex, edgeBound[edge.head]);
            ++edgeBound[edge.head];

            swapAdjacencyEdges(edge.tail, edge.tailIndex, edgeBound[edge.tail]);
            ++edgeBound[edge.tail];

            edge.alive = true;
            ++activeEdgeCount;
        }

        void swapnodes(int first, int second) noexcept {
            if (first == second) {
                return;
            }

            std::swap(nodes[first], nodes[second]);
            nodeIndex[nodes[first]] = first;
            nodeIndex[nodes[second]] = second;
        }

        void restorenode(int u) noexcept {
            swapnodes(nodeIndex[u], nodeBound);
            ++nodeBound;
        }

        void rollbackUnchecked(int checkpoint) {
            while (undoLog.size() > checkpoint) {
                const Undo undo = undoLog.back();
                undoLog.pop_back();

                if (undo.kind == UndoKind::node) {
                    restorenode(static_cast<int>(undo.id));
                } else if (undo.kind == UndoKind::edge) {
                    restoreEdge(undo.id);
                } else {
                    storedEdges[undo.id].data = std::move(oldData.back());
                    oldData.pop_back();
                }
            }
        }

    };

    template<typename Data>
    struct Graph<Data>::ReindexedSubgraph {
        Graph<Data> graph;

        // newToOldVertex[newVertex] = originalVertex
        std::vector<int> newToOldVertex;

        // oldToNewVertex[originalVertex] = newVertex, or -1
        std::vector<int> oldToNewVertex;

        // newToOldEdge[newEdgeId] = originalEdgeId
        std::vector<int> newToOldEdge;
    };
}
#endif //OBLIVIOUSROUTING_GRAPH_H