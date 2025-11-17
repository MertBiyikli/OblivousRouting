//
// Created by Mert Biyikli on 10.06.25.
//

#include "mcct_derandomized_weighted_solver.h"
#include <unordered_set>
#include <list>
#include "../../../utils/hash.h"
#include <cmath>

struct PairHash {
    std::size_t operator()(const std::pair<int,int>& p) const noexcept {
        return std::hash<int>()(p.first) ^ (std::hash<int>()(p.second) << 1);
    }
};

FRT_Tree MCCT_Solver::getBestTree(bool debug) {
    computeBestBetaAndPermutation();
    computeBestTree(debug);
    return this->tree;
}


void MCCT_Solver::computeBestTree(bool debug) {
    if(!this->graph->IsDistanceMatrixComputed()) {
        this->graph->createDistanceMatrix();
    }

    tree = FRT_Tree();
    double diameter = graph->GetDiameter();
    if(diameter == 0) {
        diameter = 1;
    }

    int i = static_cast<int>(std::ceil(std::log2(diameter) / std::log2(2.0))) + 1;
    if(debug ) {
        std::cout << "Diameter: " << diameter << ", i: " << i << std::endl;
    }

    tree.createRoot(i);
    tree.GetRoot()->addVertices(graph->getVertices());
    tree.GetRoot()->setCenter(verticesPermutation.front()); // Set the center of the root node to vertex 0
    i--;

    while(i >= 0) {
        double betaI = std::pow(2.0, i-1)*bestBeta;
        if(debug) {
            std::cout << "Processing layer " << i << " with beta: " << betaI << std::endl;
        }

        std::set<int> assigned;
        for(int& v : verticesPermutation) {
            if(assigned.find(v) != assigned.end()) {
                continue; // Skip if the vertex is already assigned
            }

            for(std::shared_ptr<FRT_Node>& node : tree.GetLayerNodes(i+1)) {
                std::shared_ptr<FRT_Node> child = std::make_shared<FRT_Node>();
                for(int vertexInCluster : node->GetVertices()) {
                    if(assigned.find(vertexInCluster) != assigned.end()) {
                        continue; // Skip if the vertex is already assigned
                    }

                    if(graph->getShortestDistance(v, vertexInCluster) < betaI) {
                        if(debug) {
                            std::cout << "Adding vertex " << vertexInCluster << " to node with center " << v << std::endl;
                        }
                        child->addVertex(vertexInCluster);
                        assigned.insert(vertexInCluster);
                    }
                }

                if(!child->GetVertices().empty()) {
                    if(debug) {
                        std::cout << "Creating child node with center " << v << " for parent node with center " << node->getCenter() << std::endl;
                    }
                    node->addChild(child);
                    child->setParent(node);
                    child->setCenter(v);

                    tree.AddNodeToLayer(i, child);
                }

            }
        }
        i--;
    }
}



void MCCT_Solver::computeBestBetaAndPermutation() {
    this->graph->createDistanceMatrix();
    this->computeBetas(debug);


    double minExpectation = std::numeric_limits<double>::max();
    std::set<int> unsettledVertices(graph->getVertices().begin(), graph->getVertices().end());

    for(double beta : betas) {
        double expectedCost = computeExpectation(beta, verticesPermutation, unsettledVertices, debug);
        if (expectedCost < minExpectation) {
            bestBeta = beta;
            minExpectation = expectedCost;
        }
    }

    while(this->verticesPermutation.size() < graph->getNumNodes()) {
        // Find the vertex with the minimum
        int bestVertex = -1;
        minExpectation = std::numeric_limits<double>::max();
        std::vector<int> candidates(this->verticesPermutation);
        for(int v : unsettledVertices) {
            if(std::find(candidates.begin(), candidates.end(), v) != candidates.end()) {
                continue; // Skip if the vertex is already in the permutation
            }

            candidates.push_back(v);
            double expectation = computeExpectation(bestBeta, candidates, unsettledVertices, debug);

            if(expectation < minExpectation) {
                minExpectation = expectation;
                bestVertex = v;
            }
            candidates.pop_back(); // Remove the vertex after checking
        }
        this->verticesPermutation.push_back(bestVertex);
        unsettledVertices.erase(bestVertex); // Remove the vertex from unsettled vertices
        this->removeDemands(bestBeta, bestVertex);
    }

}


void MCCT_Solver::computeBetas(bool debug) {
    auto& vertices = graph->getVertices();
    std::list<int> toBeAnalyzed(vertices.begin(), vertices.end());

    for(int& v : vertices) {
        toBeAnalyzed.pop_front(); // Remove vertex 0 from the list
        for(int& u : toBeAnalyzed) {
            // add debug statement to check the vertices being analyzed
            if (v == u) continue; // Skip if both vertices are the same
            // Calculate the beta value based on the distance between vertices v and u
            if (v < 0 || u < 0 || v >= vertices.size() || u >= vertices.size()) {
                throw std::out_of_range("Vertex index out of range");
            }


            double distance = graph->getShortestDistance(v, u);

            if (distance < 0) {
                throw std::runtime_error("Edge distance must be non-negative.");
            }

            double log2Dist = std::log2(distance);
            double base = std::pow(2.0, std::floor(log2Dist));
            double beta = distance / base;

            if(debug) {
                std::cout << "Distance between vertices " << v << " and " << u << ": " << distance << std::endl;
                std::cout << "Log2 distance: " << log2Dist << ", Base: " << base << ", Beta: " << beta << std::endl;
            }

            if(std::isnan(beta) || std::isinf(beta)) {
                throw std::runtime_error("Beta value is NaN or Inf, check the distance values.");
            }
            // Enforce minimum beta > 1
            if (beta <= 1.0) beta = 2.0;

            betas.insert(beta);

            // Add the beta value to the set
            if(debug)
                std::cout << "Beta for vertices " << v << " and " << u << ": " << beta << std::endl;
            // Add the beta value to the set

        }
    }
}

double MCCT_Solver::computeExpectation(double beta, std::vector<int> &verticesPermutation,
                                       std::set<int> &unsettledVertices, bool debug) {
    double result = 0.0;
    std::unordered_set<std::pair<int, int>, PairHash> settledDemands;

    if(!verticesPermutation.empty()) {
        int center = verticesPermutation.back();

        for(const auto&[u, demandMap] : this->idVertex2idVertex2demand) {
            for(const auto&[v, d_v] : demandMap) {
                if(debug) {
                    std::cout << "Checking demand from " << u << " to " << v << " with demand " << d_v << std::endl;
                }

                double firstDistance = graph->getShortestDistance(center, u);
                double secondDistance = graph->getShortestDistance(center, v);

                double bigger = std::max(firstDistance, secondDistance);
                double smaller = std::min(firstDistance, secondDistance);

                bigger /= beta;
                smaller /= beta;

                int powerBigger = static_cast<int>(std::ceil(std::log2(bigger)/std::log2(2.0))+2);
                int powerSmaller = static_cast<int>(std::ceil(std::log2(smaller)/std::log2(2.0))+2);

                if(debug)
                    std::cout << "Power bigger: " << powerBigger << ", Power smaller: " << powerSmaller << std::endl;

                if(powerBigger != powerSmaller) {
                    result += d_v * (std::pow(2.0, powerBigger+2));
                    settledDemands.insert({u, v});
                }else{
                    int level = this->demand2levelIncluded[{u, v}];
                    if (debug) {
                        std::cout << "Current level for demand from " << u << " to " << v << ": " << level << std::endl;
                    }

                    if(level == 0) {
                        level = std::numeric_limits<int>::max();
                    }
                    if(powerBigger<level){
                        this->demand2levelIncluded[{u, v}] = static_cast<int>(powerBigger);

                        if(debug) {
                            std::cout << "Updating level for demand from " << u << " to " << v << " to " << powerBigger << std::endl;
                        }
                    }
                }

            }
        }
    }

    // compute the expected cost for the unsettled vertices
    for(const auto&[u, demand_map] : this->idVertex2idVertex2demand) {
        for(const auto&[v, d_uv] : demand_map) {
            if(settledDemands.contains({u, v})) {
                if(debug) {
                    std::cout << "Demand from " << u << " to " << v << " is already settled." << std::endl;
                }
                continue; // Skip already settled demands
            }

            const int level_u_v = demand2levelIncluded[{u, v}];

            std::unordered_map<int, double> centerCosts;
            for(const auto& center : unsettledVertices)
            {
                if(!verticesPermutation.empty()
                && verticesPermutation.back() == center) {
                    continue; // Skip the last vertex in the permutation
                }

                const double& firstDistance = graph->getShortestDistance(center, u);
                const double& secondDistance = graph->getShortestDistance(center, v);

                double bigger = std::max(firstDistance, secondDistance);
                double smaller = std::min(firstDistance, secondDistance);

                bigger /= beta;
                smaller /= beta;

                double powerBigger = static_cast<int>(std::ceil(std::log2(bigger)/std::log2(2.0)))+2;
                double powerSmaller = static_cast<int>(std::ceil(std::log2(smaller)/std::log2(2.0)))+2;

                if(debug) {
                    std::cout << "Center: " << center << ", Power bigger: " << powerBigger
                              << ", Power smaller: " << powerSmaller << std::endl;
                }

                if(powerBigger != powerSmaller) {
                    if(powerSmaller == std::numeric_limits<int>::min()+2) {
                        powerSmaller = 0; // Handle the case where powerSmaller is too small
                    }

                    for(int i = powerBigger; i > powerSmaller; --i) {
                        if(level_u_v < i) {
                            centerCosts[center] += d_uv * std::pow(2.0, i+2);
                            break; // Break after adding the cost for this center
                        }
                    }
                }
            }
            for(const auto& [key, value]:centerCosts) {
                result += value/centerCosts.size();
            }

        }
    }
    return result;
}

void MCCT_Solver::removeDemands(double beta, int _v) {
    std::unordered_map<int, std::vector<int>> idVertex2vertices;
    int center = _v;
    for(auto&[u, demand_map] : this->idVertex2idVertex2demand) {
        auto& vertices = idVertex2vertices[u];
        for(auto& [v, d_uv] : demand_map) {
            double firstDistance = graph->getShortestDistance(center, u);
            double secondDistance = graph->getShortestDistance(center, v);

            double bigger = std::max(firstDistance, secondDistance);
            double smaller = std::min(firstDistance, secondDistance);

            bigger /= beta;
            smaller /= beta;

            int powerBigger = static_cast<int>(std::floor(std::log(bigger) / std::log(2))) + 2;
            int powerSmaller = static_cast<int>(std::floor(std::log(smaller) / std::log(2))) + 2;

            if (powerSmaller != powerBigger) {
                if (vertices.empty()) {
                    vertices = {v};
                }else{
                    vertices.push_back(v);
                }
            }
        }
    }

    for(auto& [u, vertices] : idVertex2vertices) {
        for(auto& v : vertices) {
            idVertex2idVertex2demand[u].erase(v);
            if(idVertex2idVertex2demand[u].empty()) {
                idVertex2idVertex2demand.erase(u);
            }
        }
    }
}
















/*

MCCTDerandomizedWeightedSolver::MCCTDerandomizedWeightedSolver() {
    this->betas = std::set<double>();
    this->verticesPermutation = std::vector<int>();
    this->idVertex2idVertex2demand = std::map<int, std::map<int, double>>();
    this->demand2levelIncluded = std::map<std::pair<int, int>, int>();
}

// Add demand from vertex v1 to vertex v2
void MCCTDerandomizedWeightedSolver::addDemand(int v1, int v2, double demand) {
    if(v1 < 0 || v2 < 0) {
        throw std::invalid_argument("Vertex IDs must be non-negative.");
    }
    if (v1 == v2) {
        throw std::invalid_argument("Cannot add demand from a vertex to itself.");
    }
    // Ensure that v1 and v2 are in the idVertex2idVertex2demand map
    if (idVertex2idVertex2demand.find(v1) == idVertex2idVertex2demand.end()) {
        idVertex2idVertex2demand[v1] = std::map<int, double>();
    }

    this->idVertex2idVertex2demand[v1][v2] = demand;
}


Tree MCCTDerandomizedWeightedSolver::getTree() const {
    return tree;
}

void MCCTDerandomizedWeightedSolver::setTree(const Tree& tree) {
    this->tree = tree;
}

std::shared_ptr<Graph> MCCTDerandomizedWeightedSolver::getGraph() const {
    return graph;
}

void MCCTDerandomizedWeightedSolver::setGraph(std::shared_ptr<Graph> graph) {
    this->graph = graph;
}


std::set<double> MCCTDerandomizedWeightedSolver::getBetas() const {
    return betas;
}

void MCCTDerandomizedWeightedSolver::setBetas(const std::set<double>& betas) {
    this->betas = betas;
}

double MCCTDerandomizedWeightedSolver::getBestB() const {
    return bestB;
}

void MCCTDerandomizedWeightedSolver::setBestB(double bestB) {
    this->bestB = bestB;
}

std::vector<int> MCCTDerandomizedWeightedSolver::getVerticesPermutation() const {
    return verticesPermutation;
}

void MCCTDerandomizedWeightedSolver::setVerticesPermutation(const std::vector<int>& verticesPermutation) {
    this->verticesPermutation = verticesPermutation;
}

Tree MCCTDerandomizedWeightedSolver::getBestTree() {
    findBestBetaAndPermutation();
    computeBestTree();
    return tree;
}


void MCCTDerandomizedWeightedSolver::findBestBetaAndPermutation() {
    // 1. Compute betas based on demands/distances
    graph->createDistanceMatrix();
    computeBetas();


    // 2. init verticesPermutation and unsettledVertices
    double minExpectation = std::numeric_limits<double>::max();
    std::set<int> unsettledVertices;
    for(const auto& i : graph->GetVertices() ) {
        unsettledVertices.insert(i);
    }


    // 3. Compute best beta based on expected costs
    for (double beta : betas) {
        double expectedCost = computeExpectation(beta, verticesPermutation, unsettledVertices);
        if (expectedCost < minExpectation) {
            bestB = beta;
            minExpectation = expectedCost;
        }
    }

    // 4. Greedily compute best permutation of the vertices
    while (verticesPermutation.size() < graph->GetVertices().size()) {
        int bestV = 0;
        minExpectation = std::numeric_limits<double>::max();
        std::vector<int> candidatePermutation = verticesPermutation;
        for (const auto& candidateVertex : unsettledVertices) {
            candidatePermutation.push_back(candidateVertex);
            double expectedCost = computeExpectation(bestB, candidatePermutation, unsettledVertices);
            if (expectedCost <= minExpectation) {
                bestV = candidateVertex;
                minExpectation = expectedCost;
            }
            candidatePermutation.pop_back();
        }
        verticesPermutation.push_back(bestV);
        unsettledVertices.erase(bestV);
        removeDemands(bestB, bestV);
    }
}


void MCCTDerandomizedWeightedSolver::computeBestTree() {
    if (graph->GetDistanceMatrix().empty()) {
        graph->createDistanceMatrix();
    }

    tree = Tree();
    double diameter = graph->GetGraphDiameter();
    if (diameter == 0) diameter = 1;
    int i = static_cast<int>(std::ceil(std::log(diameter) / std::log(2))) + 1;
    //int maxDepth = static_cast<int>(std::ceil(std::log(diameter) / std::log(2))) + 1;

    tree.CreateRoot(i);
    tree.GetRoot()->AddVertices(graph->GetVertices());
    tree.GetRoot()->SetCenter(verticesPermutation.front());

    i--;

    // iterate through the layers of the tree
    int AllLeavesConstructed = 0;
    std::unordered_set<int> assignedVerticesId;
    //while( ! tree.IsLayerWithAllSingletons(i)) {
    while (i >= 0) {

        double betaI = std::pow(2, i - 1) * bestB;
        for (auto& vertex : verticesPermutation) {

            if(assignedVerticesId.contains(vertex)) continue;
            for (auto& node : tree.GetLayerNodes(i + 1)) {

                if( node->GetVertices().size() <= 1 ) {
                    continue;
                }
                std::shared_ptr<Node> childNode = std::make_shared<Node>();

                for (auto &vertexInCluster: node->GetVertices()) {
                    if (assignedVerticesId.contains(vertexInCluster)) continue;

                    double distance = graph->GetDistanceMatrix().at(vertex).at(vertexInCluster);
                    if (distance < betaI) {
                        childNode->AddVertex(vertexInCluster);
                        assignedVerticesId.insert(vertexInCluster);
                    }
                }

                if (!childNode->GetVertices().empty()) {

                    childNode->SetParent(node);
                    childNode->SetCenter(childNode->GetVertices().front()); // TODO find a proper way of setting the center
                    // Move ownership into the tree:
                    node->AddChild(childNode);

                    if (childNode->GetVertices().size() == 1) {
                        AllLeavesConstructed++;
                    }
                    //treeArcMap_[{node->GetCenter(), childNode->GetCenter()}] = getArcBetween(node->GetCenter(), childNode->GetCenter());

                    // Now childNode == nullptr, but we already saved rawChildPtr above.
                    tree.AddNodeToLayer(i, childNode);
                }
            }

        }
        assignedVerticesId.clear();
        --i;
    }
}

void MCCTDerandomizedWeightedSolver::removeDemands(double beta, int bestV) {
    std::map<int, std::vector<int>> idVertex2vertices;
    auto& idVertex2distance = graph->GetDistanceMatrix().at(bestV);

    for (const auto& pair : idVertex2idVertex2demand) {
        int idVertex = pair.first;
        auto& idVertex2Demand = pair.second;
        std::vector<int> vertices;

        for (const auto& demandPair : idVertex2Demand) {
            int idVertex2 = demandPair.first;
            double firstDistance = idVertex2distance.at(idVertex);
            double secondDistance = idVertex2distance.at(idVertex2);
            double bigger = std::max(firstDistance, secondDistance);
            double smaller = std::min(firstDistance, secondDistance);
            bigger /= beta;
            smaller /= beta;
            int powerBigger = static_cast<int>(std::floor(std::log(bigger) / std::log(2)) + 2);
            int powerSmaller = static_cast<int>(std::floor(std::log(smaller) / std::log(2)) + 2);

            if (powerSmaller != powerBigger) {
                if (vertices.empty()) {
                    idVertex2vertices[idVertex] = {idVertex2};
                } else {
                    vertices.push_back(idVertex2);
                }
            }
        }
    }

    for (const auto& pair : idVertex2vertices) {
        int idVertex = pair.first;
        for (int idVertex2 : pair.second) {
            idVertex2idVertex2demand[idVertex].erase(idVertex2);
            if (idVertex2idVertex2demand[idVertex].empty()) {
                idVertex2idVertex2demand.erase(idVertex);
            }
        }
    }
}


double MCCTDerandomizedWeightedSolver::computeExpectation(double beta, std::vector<int>& verticesPermutation, std::set<int>& unsettledVertices, bool debug) {
    double result = 0.0;
    std::unordered_set< std::pair<int, int> > cutAndSettled;

    if (!verticesPermutation.empty()) {
        int center = verticesPermutation.back();
        auto& idVertex2distance = graph->GetDistanceMatrix()[center];

        for (const auto& pair : idVertex2idVertex2demand) {
            int idVertex = pair.first;
            auto& idVertex2Demand = pair.second;

            for (const auto& demandPair : idVertex2Demand) {
                int idVertex2 = demandPair.first;
                double demand = demandPair.second;

                double firstDistance = idVertex2distance[idVertex];// != idVertex2distance.end() ? idVertex2distance.at(idVertex) : std::numeric_limits<double>::infinity();
                double secondDistance = idVertex2distance[idVertex2];// != idVertex2distance.end() ? idVertex2distance.at(idVertex2) : std::numeric_limits<double>::infinity();
                double bigger = std::max(firstDistance, secondDistance);
                double smaller = std::min(firstDistance, secondDistance);
                bigger /= beta;
                smaller /= beta;
                int powerBigger = static_cast<int>(std::floor(std::log(bigger) / std::log(2)) + 2);
                int powerSmaller = static_cast<int>(std::floor(std::log(smaller) / std::log(2)) + 2);

                if (powerSmaller != powerBigger) {
                    if(debug) {
                        std::cout << "Center " << center << " cuts " << idVertex << " and " << idVertex2
                                  << " with cost " << std::pow(2, powerBigger + 2) * demand << " at level "
                                  << powerBigger - 1 << std::endl;
                    }
                    result += std::pow(2, powerBigger + 2) * demand;
                    cutAndSettled.insert({idVertex, idVertex2});
                }else{
                    auto demandKey = std::make_pair(idVertex, idVertex2);
                    int previousLevel = (demand2levelIncluded.count(demandKey))
                                        ? demand2levelIncluded.at(demandKey)
                                        : std::numeric_limits<int>::max();

                    if(powerBigger < previousLevel) {
                        demand2levelIncluded[{idVertex, idVertex2}] = powerBigger;
                    }
                }
            }
        }
        if(debug) {
            std::cout << "Center " << center << " has " << cutAndSettled.size() << " settled demands." << std::endl;
        }
    }

    // compute expected cost
    for(const auto& [v, v2Demand] : idVertex2idVertex2demand) {
        for(const auto& [v2, demand] : v2Demand) {

            auto demandKey = std::make_pair(v, v2);

            if(cutAndSettled.find(demandKey) != cutAndSettled.end()) {
                continue;
            }

            std::unordered_map<int, double> center2Cost;
            for(int center: unsettledVertices) {
                if(!verticesPermutation.empty() && center == verticesPermutation.back()) {
                    continue; // Skip the last vertex in the permutation
                }

                const auto& centerDistanceMap = graph->GetDistanceMatrix()[center];
                double firstDist = centerDistanceMap.at(v);
                double secondDist = centerDistanceMap.at(v2);

                double bigger = std::max(firstDist, secondDist);
                double smaller = std::min(firstDist, secondDist);

                int powerBigger = static_cast<int>(std::floor(std::log(bigger) / std::log(2))) + 2;
                int powerSmaller = static_cast<int>(std::floor(std::log(smaller) / std::log(2))) + 2;

                // special case for very small values
                if (powerSmaller == std::numeric_limits<int>::min() + 2)
                    powerSmaller = 0;

                if (powerSmaller != powerBigger) {
                    for (int i = powerBigger; i > powerSmaller; --i) {
                        if (demand2levelIncluded.find(demandKey) != demand2levelIncluded.end() && demand2levelIncluded[demandKey] < i) {
                            center2Cost[center] = std::pow(2, i + 2) * demand;
                            break;
                        }
                    }
                }
            }
            // Add averaged cost from all centers that could separate this pair
            for (const auto& [_, cost] : center2Cost) {
                result += cost / center2Cost.size();
            }
        }
    }

    return result;
}


void MCCTDerandomizedWeightedSolver::computeBetas() {
    const auto& vertices = graph->GetVertices();
    const auto& distanceMatrix = graph->GetDistanceMatrix();

    for (size_t i = 0; i < vertices.size(); ++i) {
        int u = vertices[i];
        for (size_t j = i + 1; j < vertices.size(); ++j) {
            int v = vertices[j];

            double distance = distanceMatrix[u][v];
            if (distance <= 0.0 || std::isnan(distance)) continue;

            double log2Dist = std::log2(distance);
            double base = std::pow(2.0, std::floor(log2Dist));
            double beta = distance / base;

            // Enforce minimum beta > 1
            if (beta <= 1.0) beta = 2.0;

            betas.insert(beta);
        }
    }
}



void MCCTDerandomizedWeightedSolver::reset() {
    betas.clear();
    verticesPermutation.clear();
    idVertex2idVertex2demand.clear();
    tree.clear();
    graph = nullptr;
}

*/