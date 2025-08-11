
#include "raecke_tree_decomp.h"
#include "mcct/mcct_derandomized_weighted_solver.h"

#include <queue>
#include <unordered_set>
#include <limits>
#include <cmath>
#include <iostream>

/*
void RackeObliviousRoutingSolver::createTreesAndLambda() {
    // If the solver object doesn’t exist yet, build it and give it our graph_
    if (!mcct_) {
        mcct_ = std::make_unique<MCCTDerandomizedWeightedSolver>();
    }
    if (!mcct_->getGraph()) {
        // If the graph is not set, set it now
        mcct_->setGraph(graph_);
    }

    double lambdaSum = 0.0;
    int treeIndex = 0;

    // Keep building until sum of lambdas reaches 1:
    while (lambdaSum < 1.0) {
        // 1) Compute a new tree from the MCCT solver
        auto copyGraph = std::make_shared<Graph>(*graph_);
        graphs_.push_back(copyGraph);

        if (!mcct_) {
            mcct_ = std::make_unique<MCCTDerandomizedWeightedSolver>();
        }

        auto t = getTree(copyGraph);
        //mcct_.reset();
        trees_.push_back(t);

        // 3) Build demands in that copy so we can compute r‐loads later
        // setRequirements(copyGraph);

        // 4) Compute R‐loads for this treeIndex
        computeRLoads(treeIndex, t, copyGraph);

        // 5) Compute max R‐load
        double rmax = getMaxRLoad(treeIndex, t);

        // 6) Decide how much λᵢ to assign
        // TODO: why should the second argument be 1.0/rmax even though in the pseudocode this is rmax only
        double l = std::min(1-lambdaSum, 1.0 / rmax);
        lambdas_.push_back(l);
        lambdaSum += l;

        treeIndex++;
    }
}


void RackeObliviousRoutingSolver::setRequirements(const std::shared_ptr<Graph>& g) {

    for (const auto& uID : g->GetVertices()) {
        for (auto& arcPtr : g->neighbors(uID)) {
            int vID = arcPtr->target;
            double cap = arcPtr->capacity;

            mcct_->addDemand(uID, vID, cap);
        }
    }
}


// TODO: compute the rloads as defined in the paper by summing up the capacity of each tree edge that lies on the unique
//  path in the tree. So basically if you want to compute the load for the edge (u, v) in the graph,
//  you have to traverse from one leaf node of the tree to the other leaf node, and sum up the capacities of all tree edges which basically
//  "builds" up your cut from each "clustering"
void RackeObliviousRoutingSolver::computeRLoads(int treeIndex,
                                                const std::shared_ptr<Tree>& t,
                                                const std::shared_ptr<Graph>& copyGraph)
{

    idTree2arc2rLoad_.emplace_back();
    auto& arc2rLoad = idTree2arc2rLoad_.back();

    std::queue<std::shared_ptr<Node>> bfsQueue;
    auto rootPtr = t->GetRoot();
    for (auto& childPtr : rootPtr->GetChildren()) {
        bfsQueue.push(childPtr);
    }

    // While queue not empty, pop front:
    while (!bfsQueue.empty()) {
        auto nodePtr = bfsQueue.front();
        bfsQueue.pop();

        // Parent is a weak_ptr<Node>, so lock it:
        auto parentPtr = nodePtr->GetParent();
        if (parentPtr) {

            int centerParent = parentPtr->GetCenter();
            int currentCenter = nodePtr->GetCenter();

            // If the node are the same, we do not need to store the cut sum
            if (parentPtr->GetVertices().size() == nodePtr->GetVertices().size()) {
                // Enqueue grandchildren:
                for (auto &gc: nodePtr->GetChildren()) {
                    bfsQueue.push(gc);
                }
                continue;
            }



            std::vector<int> clusterVerts = nodePtr->GetVertices();
            std::unordered_set<int> remaining(copyGraph->GetVertices().begin(), copyGraph->GetVertices().end());
            // Remove the cluster vertices from the remaining set:
            for (int vID: clusterVerts) {
                remaining.erase(vID);
            }


            // Accumulate the cut‐sum:
            double cutSum = 0.0;
            for (int uID: clusterVerts) {
                // Suppose copyGraph->getOutgoingArcs(int) returns vector<Arc*> from that uID

                for (auto &arcPtr: copyGraph->neighbors(uID)) {
                    // Only count arcs “emanating” from uID → wID
                    int wID = arcPtr->target;
                    if (remaining.contains(wID)) {
                        cutSum += arcPtr->capacity;
                    }
                }
            }

            // we have to ensure that whenever the nodes are different, w.r.t. to the contained vertices, the centers cannot be the same, since
            // otherwise we could not compute the shortest path

            if (currentCenter == centerParent) {
                // take another center vertex from the parent node
                std::set<int> otherCenters;
                for(const int& v : parentPtr->GetVertices()) {
                    if (!nodePtr->ContainsVertex(v)) {
                        centerParent = v;
                        break;
                    }
                }
            }

            auto shortestPath = copyGraph->getShortestPath(
                    centerParent,
                    currentCenter
            );


            for(auto& edge : shortestPath) {
                // Normalize the edge distance based on the cutSum
                double rload = arc2rLoad[edge];
                rload += cutSum / edge->capacity; // Normalize by capacity
                arc2rLoad[edge] = rload;
                auto rev_arc = copyGraph->GetReverseEdge(edge);
                arc2rLoad[rev_arc] = rload; // Also update the reverse arc
            }


            std::cout << "Cut sum for nodes ";
            for (const int& node : nodePtr->GetVertices()) {
                std::cout << node << " ";
            }

            std::cout << "with parent set:";
            for(const int& parent : parentPtr->GetVertices()) {
                std::cout << " " << parent;
            }
            std::cout << " is: " << cutSum << std::endl;

            // debug the edge weights
            std::cout << "Edge weights for tree index " << treeIndex << ": ";
            for (const auto& [edge, rLoad] : arc2rLoad) {
                 std::cout << "Edge (" << (edge->source) << " →  " << edge->target << ") has r-load: " << rLoad << ", ";
            }



            // Enqueue grandchildren:
            for (auto &gc: nodePtr->GetChildren()) {
                bfsQueue.push(gc);
            }
        }
    }
}


double RackeObliviousRoutingSolver::getMaxRLoad(int treeIndex,
                                                    const std::shared_ptr<Tree>& t)
{
    // Grab the r‐load map for this tree
    const auto& arc2rLoad = idTree2arc2rLoad_.at(treeIndex);
    double maxRatio = 0.0;
    for(const auto& [arc, rLoad] : arc2rLoad) {
        if (rLoad > maxRatio) {
            maxRatio = rLoad;
        }
    }

    return maxRatio;
}

RackeObliviousRoutingSolver::~RackeObliviousRoutingSolver() {

}

double RackeObliviousRoutingSolver::getRloadAllEdges(const std::shared_ptr<Graph>& g) const {
    double totalRLoadAllEdges = 0.0;
    for (const auto& v : g->GetVertices()) {
        for (const auto& u : g->neighbors(v)) {
            double totalPerEdge = 0.0;

            // Sum over all trees i:  rLoad(i,arc) * lambda[i]
            for (size_t i = 0; i < trees_.size(); ++i) {
                const auto& arcLoadMap = idTree2arc2rLoad_[i];

                auto it = arcLoadMap.find(u);
                double rLoad = (it == arcLoadMap.end() ? 0.0 : it->second);
                totalPerEdge += rLoad * lambdas_[i];
            }

            // Add e^(totalPerEdge) to the grand total
            totalRLoadAllEdges += std::exp(totalPerEdge);
        }
    }
    return totalRLoadAllEdges;
}

std::vector<std::unordered_map<std::shared_ptr<Edge>, double>> RackeObliviousRoutingSolver::GetEdge2Load() const {
    return idTree2arc2rLoad_;
}

void RackeObliviousRoutingSolver::solve(const Graph &graph) {
    std::shared_ptr<Graph> g(std::make_shared<Graph>(graph));
    this->setGraph(g);
    this->createTreesAndLambda();
}


std::vector<std::shared_ptr<Graph>>& RackeObliviousRoutingSolver::getGraphs() {
    return graphs_;
}

void RackeObliviousRoutingSolver::setGraph(const std::shared_ptr<Graph>& g) {
    graph_ = g;
}

std::vector<double> RackeObliviousRoutingSolver::getLambdas() const {
    return lambdas_;
}

std::vector<std::shared_ptr<Tree>> RackeObliviousRoutingSolver::getTrees() const{
    return trees_;
}

std::shared_ptr<Tree> RackeObliviousRoutingSolver::getTree(std::shared_ptr<Graph>& g) {
    mcct_->setGraph(g);
    setRequirements(g);
    computeNewDistances(g);

    return std::make_shared<Tree>(mcct_->getBestTree());
}

std::unordered_map<std::shared_ptr<Edge>, double> RackeObliviousRoutingSolver::computeNewDistances(std::shared_ptr<Graph>& g) {
    double totalRLoadsAllEdges = getRloadAllEdges(g);
    std::unordered_map<std::shared_ptr<Edge>, double> arc2scaledDist;

    for(const auto& v : g->GetVertices()) {
        for (const auto& edge : g->neighbors(v)) {
            double totalrLoad = 0.0;
            for (size_t i = 0; i < trees_.size(); ++i) {


                double rLoad = (idTree2arc2rLoad_[i].find(edge) != idTree2arc2rLoad_[i].end() ? idTree2arc2rLoad_[i].at(edge) : 0.0);

                totalrLoad += rLoad * lambdas_[i];
            }
            std::cout << "Total rload: " << totalrLoad << std::endl;
            double num         = std::exp(totalrLoad) / edge->capacity;
            double newDistance = num / totalRLoadsAllEdges;

            if(isinf(newDistance)
            || isinf(num)) {
                throw std::runtime_error("Infinity encountered in newDistance or Exponential calculation for edge: " + std::to_string(v) + " → " + std::to_string(edge->target));
            }

            if(isnan(newDistance)) {
                throw std::runtime_error("NaN encountered in newDistance calculation for edge: " + std::to_string(v) + " → " + std::to_string(edge->target));
            }

            arc2scaledDist[edge] = newDistance;

            //DEbug for printing
            std::cout << "Edge (" << v << " → " << edge->target << ") has new distance: "
                      << std::fixed << std::setprecision(6) << newDistance << std::endl;
            std::cout << "Rloads for edge (" << v << " → " << edge->target << "): ";
            for (size_t i = 0; i < trees_.size(); ++i) {
                double rLoad = (idTree2arc2rLoad_[i].find(edge) != idTree2arc2rLoad_[i].end() ?
                                idTree2arc2rLoad_[i].at(edge) : 0.0);
                std::cout << rLoad * lambdas_[i] << " ";
            }
        }
    }

    // Normalize distances
    normalizeDistance(g, arc2scaledDist);
    return arc2scaledDist;
}

void RackeObliviousRoutingSolver::normalizeDistance(std::shared_ptr<Graph>& _g, std::unordered_map<std::shared_ptr<Edge>, double>& arc2scaledDist) {
    double minScaled = std::numeric_limits<double>::infinity();
    for(const auto& v : _g->GetVertices()) {
        for (const auto& arc : _g->neighbors(v)) {
            double scaledDist = arc2scaledDist[arc];
            if (scaledDist < minScaled) {
                minScaled = scaledDist;
            }
        }
    }

    for(const auto& v : _g->GetVertices()) {
        for (auto& arc : _g->neighbors(v)) {
            double arc2scaledDistValue = arc2scaledDist[arc];
            double newDistance = arc2scaledDistValue/minScaled;

            if(isnan( newDistance) ) {
                bool found = false;
                throw std::runtime_error("Found NaN in newDistance for arc: " + std::to_string(v) + " → " + std::to_string(arc->target));
            }
            // TODO: recheck if this is really necessary...
            if(newDistance < 1) {
                newDistance = 1.0;
            }
            arc->distance = newDistance;
        }
    }
}

void RackeObliviousRoutingSolver::printObliviousRoutingTable() const {
    // Vertex set
    const auto& verts = graph_->GetVertices();

    // For each directed arc (u -> v) in the original graph
    for (int u : verts) {
        for (const auto& edge : graph_->neighbors(u)) {
            int v = edge->target;
            std::cout << "Arc (" << u << " → " << v << "): ";
            bool firstEntry = true;

            // For each commodity (s -> t)
            for (int s : verts) {
                for (int t : verts) {
                    if (s == t) continue;
                    // Expected use = sum_i λ_i * 1{(u->v) lies on path in tree i}
                    double expectedUse = 0.0;
                    for (size_t i = 0; i < trees_.size(); ++i) {
                        //auto treePath = trees_[i]->getPath(s, t);
                        auto path = graphs_[i]->getShortestPath(s, t);
                        bool onPath = std::any_of(
                                path.begin(), path.end(),
                                [&](const auto& ePtr){ return ePtr == edge; }
                        );
                        if (onPath) expectedUse += lambdas_[i];
                    }

                    if (expectedUse > 0.0) {
                        if (!firstEntry) std::cout << ", ";
                        std::cout << "[" << s << "→" << t << "]="
                                  << std::fixed << std::setprecision(6)
                                  << expectedUse;
                        firstEntry = false;
                    }
                }
            }
            std::cout << "\n";
        }
    }
}
*/




Graph RaeckeFRT::getGraph() const {
    return m_graph;
}

void RaeckeFRT::setGraph(const Graph &g) {
    m_graph = g;
}

void RaeckeFRT::run() {
    auto p = std::make_shared<Graph>(m_graph);
    this->m_mcct.setGraph(p);

    m_lambdaSum = 0.0;
    int treeIndex = 0;

    while(m_lambdaSum < 1) {
        if(debug)
            std::cout << "Computing tree " << treeIndex << std::endl;
        m_lambdaSum += iterate(treeIndex);
        treeIndex++;
    }
}

double RaeckeFRT::iterate(int treeIndex) {
    auto copyGraph = std::make_shared<Graph>(m_graph);

    FRT_Tree t = getTree(copyGraph);
    m_trees.push_back(t);
    m_mcct.reset();

    computeRLoads(treeIndex, t, copyGraph);
    double l = getMaxRload(treeIndex, t);
    double delta = std::min(1/l, 1-m_lambdaSum);
    m_lambdas.push_back(delta);
    m_lambdas.push_back(delta);
    m_graphs.push_back(*copyGraph);
    return delta;
}

FRT_Tree RaeckeFRT::getTree(std::shared_ptr<Graph> &g) {
    m_mcct.setGraph(g);
    setRequirements(g);
    computeNewDistances(g);
    return m_mcct.getBestTree();
}


void RaeckeFRT::computeRLoads(int treeIndex, FRT_Tree &_t, std::shared_ptr<Graph> &copyGraph) {
    std::queue<std::shared_ptr<FRT_Node>> bfsQueue;
    for(auto& node : _t.GetRoot()->getChildren()) {
        bfsQueue.push(node);
    }

    if(m_idTree2edge2rload.size() <= treeIndex) {
        m_idTree2edge2rload.resize(treeIndex + 1);
    }
    m_idTree2edge2rload[treeIndex] = std::unordered_map<std::pair<int, int>, double>();
    auto& edge2Load = m_idTree2edge2rload[treeIndex];

    while(!bfsQueue.empty()) {
        auto node = bfsQueue.front();
        bfsQueue.pop();
        auto parent = node->getParent();
        if(parent->getCenter() == node->getCenter()) {
            for(auto& child : node->getChildren()) {
                bfsQueue.push(child);
            }
            continue;
        }

        std::vector<int> nodeVertices;
        for (int v: node->GetVertices()) {
            nodeVertices.push_back(v);
        }
        std::vector<int> remaining(copyGraph->getVertices());
        for(int& v : nodeVertices) {
            auto it = std::find(remaining.begin(), remaining.end(), v);
            if(it != remaining.end()) {
                remaining.erase(it);
            }
        }

        double cut = 0.0;
        for(int u : nodeVertices) {
            for(const auto& v : copyGraph->neighbors(u)) {

                if(std::find(remaining.begin(), remaining.end(), v) != remaining.end()) {
                    cut += copyGraph->getEdgeCapacity(u, v);
                }
            }
        }

        int centerParent = parent->getCenter();
        int currentCenter = node->getCenter();

        auto path = copyGraph->getShortestPath(centerParent, currentCenter);
        for(int i = 0; i<path.size()-1; i++) {
            int u = path[i];
            int v = path[i+1];
            // Normalize the edge distance based on the cut
            double rLoad = edge2Load[{u, v}];
            rLoad += cut / copyGraph->getEdgeCapacity(u, v); // Normalize by capacity
            edge2Load[{u, v}] = rLoad;
            edge2Load[{v,u}] = rLoad; // Also update the reverse arc
        }
        for(auto& child : node->getChildren()) {
            bfsQueue.push(child);
        }
    }
}


double RaeckeFRT::getMaxRload(int treeIndex, FRT_Tree &_t) {
    // Grab the r‐load map for this tree
    const auto& edge2rLoad = m_idTree2edge2rload[treeIndex];
    double maxRatio = 0.0;
    for(const auto& [edge, rLoad] : edge2rLoad) {
        if (rLoad > maxRatio) {
            maxRatio = rLoad;
        }
    }

    return maxRatio;
}



void RaeckeFRT::setRequirements(const std::shared_ptr<Graph> &g) {
    for (const auto& u : g->getVertices()) {
        for (const auto& v : g->neighbors(u)) {
            double cap = g->getEdgeCapacity(u, v);
            m_mcct.addDemand(u, v, cap);
        }
    }
}

void RaeckeFRT::computeNewDistances(std::shared_ptr<Graph> &g) {
    double totalRLoadsAllEdges = this->getRloadAllEdges(g);
    std::unordered_map<std::pair<int, int>, double> edge2scaledDist;

    for(int u = 0; u<g->getNumNodes(); u++) {
        for (const auto& v : g->neighbors(u)) {
            double totalrLoad = 0.0;
            for (size_t i = 0; i < m_trees.size(); ++i) {
                auto it = this->m_idTree2edge2rload[i].find({u, v});
                double rLoad = (it != this->m_idTree2edge2rload[i].end() ? it->second : 0.0);
                totalrLoad += rLoad * m_lambdas[i];
            }
            double num         = std::exp(totalrLoad) / g->getEdgeCapacity(u, v);
            double newDistance = num / totalRLoadsAllEdges;

            if(isinf(newDistance) || isinf(num)) {
                throw std::runtime_error("Infinity encountered in newDistance or Exponential calculation for edge: " + std::to_string(u) + " → " + std::to_string(v));
            }

            if(isnan(newDistance)) {
                throw std::runtime_error("NaN encountered in newDistance calculation for edge: " + std::to_string(u) + " → " + std::to_string(v));
            }

            edge2scaledDist[{u, v}] = newDistance;
        }
    }

    // Normalize distances
    normalizeDistance(g, edge2scaledDist);
}


void RaeckeFRT::normalizeDistance(std::shared_ptr<Graph> &_g, std::unordered_map<std::pair<int, int>, double> &edge2scaledDist) {
    double minDistance = std::numeric_limits<double>::infinity();
    for(int u = 0; u < _g->getNumNodes(); u++) {
        for (const auto& v : _g->neighbors(u)) {
            double scaledDist = edge2scaledDist[{u, v}];
            if (scaledDist < minDistance) {
                minDistance = scaledDist;
            }
        }
    }

    for(int u = 0; u < _g->getNumNodes(); u++) {
        for (const auto &v: _g->neighbors(u)) {
            double arc2scaledDistValue = edge2scaledDist[{u, v}];
            double newDistance = arc2scaledDistValue / minDistance;

            if (isnan(newDistance)) {
                throw std::runtime_error(
                        "NaN encountered in newDistance for edge: " + std::to_string(u) + " → " + std::to_string(v));
            }
            if (newDistance < 1) {
                newDistance = 1.0; // Ensure minimum distance is 1}
            }
            _g->updateEdgeDistance(u, v, newDistance);
        }
    }
}


double RaeckeFRT::getRloadAllEdges(const std::shared_ptr<Graph>& g) const{
    double totalRLoadAllEdges = 0.0;
    for(int u = 0; u < g->getNumNodes(); u++) {
        for (const auto& v : g->neighbors(u)) {
            double totalPerEdge = 0.0;

            std::pair<int, int> edge = {u, v};
            // Sum over all trees i:  rLoad(i,arc) * lambda[i]
            for (size_t i = 0; i < m_trees.size(); ++i) {
                auto it = this->m_idTree2edge2rload[i].find(edge);
                if(it != this->m_idTree2edge2rload[i].end()) {
                    totalPerEdge += it->second * m_lambdas[i];
                } else {
                    // If the edge is not found, we assume the r-load is 0
                    totalPerEdge += 0.0;
                }
            }

            // Add e^(totalPerEdge) to the grand total
            totalRLoadAllEdges += std::exp(totalPerEdge);
        }
    }
    return totalRLoadAllEdges;
}
