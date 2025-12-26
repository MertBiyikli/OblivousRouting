//
// Created by Mert Biyikli on 11.12.25.
//

#include "efficient_frt.h"

#include <list>
#include <cmath>
#include <absl/strings/str_format.h>

void EfficientFRT::updateEdgeDistances(const std::vector<double> &distances) {
    for (int e = 0; e < graph.getNumEdges(); ++e) {
        graph.updateEdgeDistance(e, distances[e]);
    }
}

std::shared_ptr<EfficientFRTTreeNode> EfficientFRT::getTree() {
    for (const auto& u : graph.getVertices()) {
        for (const auto& v : graph.neighbors(u)) {
            double cap = graph.getEdgeCapacity(u, v);
            this->addDemands(u, v, cap);
        }
    }
    computeBestBetaAndPermutation();
   // auto tree = getBestTree();

    auto bottom_up = getTreeBottomUp();

    return bottom_up;
}


void EfficientFRT::cleanUpTree(std::shared_ptr<EfficientFRTTreeNode>& node) {
    if (!node) return;

    if ( node->members.size() == 1) {
        node->center = node->members.front();
    }

    // first clean children
    for (auto& ch : node->children) cleanUpTree(ch);

    // then rebuild children list safely
    std::vector<std::shared_ptr<EfficientFRTTreeNode>> new_children;
    new_children.reserve(node->children.size());

    for (auto& child : node->children) {
        if (!child) continue;

        if (same_members(node->members, child->members)) {
            // absorb child: adopt grandchildren
            for (auto& gc : child->children) {
                if (!gc) continue;
                gc->parent = node;
                new_children.push_back(gc);
            }
            child->children.clear();
            // child gets dropped
        } else {
            new_children.push_back(child);
        }
    }

    node->children.swap(new_children);
}

void EfficientFRT::computeBetas() {
    const auto& nodes = graph.getVertices();
    betas.clear();
    std::list<int> toBeAnalyzed(nodes.begin(), nodes.end());

    for (auto& v : nodes) {
        toBeAnalyzed.pop_front();
        for (auto& u : toBeAnalyzed) {
            if (v == u) continue; // Skip if both vertices are the same
            if (v < 0 || u < 0 || v >= nodes.size() || u >= nodes.size()) {
                throw std::out_of_range("Vertex index out of range");
            }

            double dist = graph.getShortestPathDistance(v, u);
            if (dist < 0.0) {
                throw std::runtime_error("Negative distance encountered between vertices " + std::to_string(v) + " and " + std::to_string(u));
            }

            double log2Dist = std::log2(dist);
            double base = std::pow(2.0, std::floor(log2Dist));
            double beta = dist / base;

            if(std::isnan(beta) || std::isinf(beta)) {
                throw std::runtime_error("Beta value is NaN or Inf, check the distance values.");
            }
            // Enforce minimum beta > 1
            if (beta <= 1.0) beta = 2.0;

            betas.insert(beta);
        }
    }
}


void EfficientFRT::computeBestBetaAndPermutation() {
    // Implementation of best beta and permutation computation
    // This is a placeholder for the actual implementation
    computeBetas();
    double minExp = std::numeric_limits<double>::infinity();
    std::unordered_set<int> allVertices(graph.getVertices().begin(), graph.getVertices().end());
    for(double beta : betas) {
        double expCost = computeExpectation(beta, allVertices, verticesPermutation);
        if (expCost < minExp) {
            bestBeta = beta;
            minExp = expCost;
        }
    }




    while(verticesPermutation.size() <
        graph.getNumNodes()) {
        // Find the vertex with the minimum
        int bestVertex = -1;
        minExp = std::numeric_limits<double>::max();
        std::vector<int> candidates(this->verticesPermutation);
        for(int v : allVertices) {
            if(std::find(candidates.begin(), candidates.end(), v) != candidates.end()) {
                continue; // Skip if the vertex is already in the permutation
            }

            candidates.push_back(v);
            double expected = computeExpectation(bestBeta, allVertices, candidates);

            if(expected < minExp) {
                minExp = expected;
                bestVertex = v;
            }
            candidates.pop_back(); // Remove the vertex after checking
        }
        this->verticesPermutation.push_back(bestVertex);
        allVertices.erase(bestVertex); // Remove the vertex from unsettled vertices
        this->removeDemands(bestBeta, bestVertex);
    }

}

std::shared_ptr<EfficientFRTTreeNode> EfficientFRT::getBestTree() {
    double diameter = graph.GetDiameter();
    if (diameter == 0.0) {
        diameter = 1;
    }

    int i = static_cast<int>(std::ceil(std::log2(diameter) / std::log2(2.0))) + 1;
    // initialize level2nodes
    level2nodes.clear();
    level2nodes.resize(i + 1);

    // Create root node
    auto root = std::make_shared<EfficientFRTTreeNode>();
    root->id = 0;
    root->members.reserve(graph.getNumNodes());
    for (int v = 0; v < graph.getNumNodes(); ++v) {
        root->members.push_back(v);
    }
    root->center = verticesPermutation[0];
    level2nodes[i].push_back(root);

    i--;


    while (i>=0) {
        double cur_beta = bestBeta * std::pow(2.0, i-1);

        std::unordered_set<int> assigned;
        for (int& v : verticesPermutation) {
            if (assigned.find(v) != assigned.end()) {
                continue;
            }

            for (auto& node : level2nodes[i+1]) {
                std::shared_ptr<EfficientFRTTreeNode> childNode = std::make_shared<EfficientFRTTreeNode>();
                for (int vertex : node->members) {
                    if (assigned.find(vertex) != assigned.end()) {
                        continue;
                    }
                    double dist = graph.getShortestPathDistance(v, vertex);

                    if (dist < cur_beta) {
                        childNode->members.push_back(vertex);
                        assigned.insert(vertex);
                    }
                }

                if (!childNode->members.empty()
                    && !(childNode->getMembers().size() == 1 && node->getMembers().size() == 1)) {

/*
                    // TODO: this needs to be fixed..., but for now
                    if (std::find(childNode->members.begin(), childNode->members.end(), v) == childNode->members.end()) {
                        childNode->center = childNode->members[0];
                    } else {
                        childNode->center = v;
                    }

                    if (node->members.size() == childNode->members.size()) {
                        childNode->center = node->center;
                    }*/
                    childNode->center = v;
                    childNode->id = static_cast<int>(level2nodes[i].size());
                    childNode->parent = node;
                    node->children.push_back(childNode);
                    level2nodes[i].push_back(childNode);
                }
            }
        }
        i--;
    }
    return root;
}


std::shared_ptr<EfficientFRTTreeNode> EfficientFRT::getTreeBottomUp() {
    const int n = graph.getNumNodes();
    if ((int)verticesPermutation.size() != n) {
        throw std::runtime_error("Vertices permutation size != n");
    }

    double diameter = graph.GetDiameter();
    if (diameter <= 0.0) diameter = 1.0;
    const int L = (int)std::ceil(std::log2(diameter)/std::log2(2));

    // Δ_i = β · 2^i   (canonical)
    auto Delta = [&](int i) -> double {
        return bestBeta * (double)(1ull << (i));
    };

    // Root cluster = V
    auto root = std::make_shared<EfficientFRTTreeNode>();
    root->center = verticesPermutation.front(); // arbitrary valid center
    root->id = -1;
    root->members.clear();
    root->members.reserve((size_t)n);
    for (int v = 0; v < n; ++v) root->members.push_back(verticesPermutation[v]);

    // Current frontier clusters at level i+1
    std::vector<std::shared_ptr<EfficientFRTTreeNode>> frontier;
    frontier.push_back(root);

    // Optional correctness checks
    const bool VERIFY_COVERAGE = true;
    const bool VERIFY_RADIUS   = false; // expensive

    auto fail = [&](const std::string& msg) {
        throw std::runtime_error("FRT verify failed: " + msg);
    };

    // Build children level by level (top-down ensures nestedness)
    for (int i = L; i >= 0; --i) {
        const double R = Delta(i);

        std::vector<std::shared_ptr<EfficientFRTTreeNode>> nextFrontier;
        nextFrontier.reserve(frontier.size() * 2);

        for (auto& parent : frontier) {
            if (!parent || parent->members.empty()) continue;

            // If already singleton, keep it as leaf.
            if (parent->members.size() == 1) {
                nextFrontier.push_back(parent);
                continue;
            }

            // Mark which vertices of this parent are still unassigned at this level
            std::unordered_set<int> unassigned;
            unassigned.reserve(parent->members.size() * 2);
            for (int v : parent->members) unassigned.insert(v);

            // Build child clusters inside this parent, scanning centers in permutation order
            parent->children.clear();

            for (int c : verticesPermutation) {
                if (unassigned.empty()) break;
                if (!unassigned.contains(c)) continue; // center must be in this parent cluster

                auto child = std::make_shared<EfficientFRTTreeNode>();
                child->center = c;
                child->parent = parent;
                child->members.clear();

                // Ball intersection: { v in unassigned ∩ parent : d(c,v) < R }
                // iterate over parent's members for locality
                for (int v : parent->members) {
                    if (!unassigned.contains(v)) continue;
                    if (graph.getShortestPathDistance(c, v) < R) {
                        child->members.push_back(v);
                    }
                }

                if (child->members.empty()) continue;

                // Remove assigned vertices from unassigned
                for (int v : child->members) unassigned.erase(v);

                parent->children.push_back(child);
                nextFrontier.push_back(child);

            }

            // Safety: every vertex in parent must be assigned (should always hold)
            if (!unassigned.empty()) {
                // This can only happen if distances behave strangely or permutation missing nodes
                // Fallback: create singleton children for leftover vertices to keep tree valid
                for (int v : unassigned) {
                    auto child = std::make_shared<EfficientFRTTreeNode>();
                    child->center = v;
                    child->id = -1;
                    child->parent = parent;
                    child->members = {v};
                    parent->children.push_back(child);
                    nextFrontier.push_back(child);
                }
                unassigned.clear();
            }

            // Coverage sanity inside this parent
            if (VERIFY_COVERAGE) {
                std::vector<char> seen(n, 0);
                size_t count = 0;
                for (auto& ch : parent->children) {
                    for (int v : ch->members) {
                        if (seen[v]) fail("duplicate vertex in children of same parent at level i=" + std::to_string(i));
                        seen[v] = 1;
                        count++;
                    }
                }
                if (count != parent->members.size()) {
                    fail("children do not cover parent.members at level i=" + std::to_string(i));
                }
            }
        }

        frontier.swap(nextFrontier);
    }


    // Ensure singleton leaves
    std::vector<std::shared_ptr<EfficientFRTTreeNode>> stack;
    stack.push_back(root);

    while (!stack.empty()) {
        auto node = stack.back();
        stack.pop_back();
        if (!node) continue;

        if (node->children.empty()
            && node->members.size() > 1) {
            // Split into singleton children
            for (int v : node->members) {
                auto leaf = std::make_shared<EfficientFRTTreeNode>();
                leaf->center = v;
                leaf->members = {v};
                leaf->parent = node;
                node->children.push_back(leaf);
            }
        } else {
            for (auto& ch : node->children) stack.push_back(ch);
        }
    }


    return root;
}




void EfficientFRT::computeCurrentPartition(std::vector<int>& centers, std::vector<int>& clusters, const double scale) {
    clusters = std::vector<int>(graph.getNumNodes(), -1);
    std::unordered_set<int> assigned;
    for (const auto& v : verticesPermutation) {
        if (assigned.find(v) != assigned.end()) {
            continue;
        }

        // v is the new center for the next cluster
        centers.push_back(v);
        assigned.insert(v);

        for (const auto& u : verticesPermutation) {
            if (assigned.find(u) != assigned.end()) {
                continue;
            }

            double dist = graph.getShortestPathDistance(v, u);
            if (dist < scale) {
                assigned.insert(u);
                clusters[u] = v;
            }
        }
    }

    // also set the cluster for centers
    for (const auto& center : centers) {
        clusters[center] = center;
    }
}


double EfficientFRT::computeExpectation(double beta, std::unordered_set<int> &allVertices, const std::vector<int> &currentPermutation) {
     double result = 0.0;
    std::unordered_set<std::pair<int, int>, PairHash> settledDemands;

    if(!currentPermutation.empty()) {
        int center = currentPermutation.back();

        for(const auto&[u, demandMap] : this->idVertex2idVertex2demand) {
            for(const auto&[v, d_v] : demandMap) {



                double firstDistance = graph.getShortestPathDistance(center, u);
                double secondDistance = graph.getShortestPathDistance(center, v);

                double bigger = std::max(firstDistance, secondDistance);
                double smaller = std::min(firstDistance, secondDistance);

                bigger /= beta;
                smaller /= beta;

                int powerBigger = static_cast<int>(std::ceil(std::log2(bigger)/std::log2(2.0))+2);
                int powerSmaller = static_cast<int>(std::ceil(std::log2(smaller)/std::log2(2.0))+2);

                if(powerBigger != powerSmaller) {
                    result += d_v * (std::pow(2.0, powerBigger+2));
                    settledDemands.insert({u, v});
                }else{
                    int level = this->demand2levelIncluded[{u, v}];


                    if(level == 0) {
                        level = std::numeric_limits<int>::max();
                    }
                    if(powerBigger<level){
                        this->demand2levelIncluded[{u, v}] = static_cast<int>(powerBigger);

                    }
                }

            }
        }
    }

    // compute the expected cost for the unsettled vertices
    for(const auto&[u, demand_map] : this->idVertex2idVertex2demand) {
        for(const auto&[v, d_uv] : demand_map) {
            if(settledDemands.contains({u, v})) {
                continue; // Skip already settled demands
            }

            const int level_u_v = demand2levelIncluded[{u, v}];

            std::unordered_map<int, double> centerCosts;
            for(const auto& center : allVertices)
            {
                if(!verticesPermutation.empty()
                && verticesPermutation.back() == center) {
                    continue; // Skip the last vertex in the permutation
                }

                const double& firstDistance = graph.getShortestPathDistance(center, u);
                const double& secondDistance = graph.getShortestPathDistance(center, v);

                double bigger = std::max(firstDistance, secondDistance);
                double smaller = std::min(firstDistance, secondDistance);

                bigger /= beta;
                smaller /= beta;

                double powerBigger = static_cast<int>(std::ceil(std::log2(bigger)/std::log2(2.0)))+2;
                double powerSmaller = static_cast<int>(std::ceil(std::log2(smaller)/std::log2(2.0)))+2;


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

void EfficientFRT::addDemands(int u, int v, double demand) {
    if (u < 0 || v < 0 || u >= graph.getNumNodes() || v >= graph.getNumNodes()) {
        throw std::out_of_range("Vertex index out of range");
    }
    if (u == v) {
        throw std::invalid_argument("Cannot add demand from a vertex to itself.");
    }
    // Ensure that u and v are in the idVertex2idVertex2demand map
    auto& it = idVertex2idVertex2demand[u];
    if (it.empty()) {
        idVertex2idVertex2demand[u] = std::unordered_map<int, double>();
    }
    it[v] = demand;
/*
    auto& v_vec = vertex2vertex[u];
    // Check if v already exists in the adjacency list
    auto it_v = std::find(v_vec.begin(), v_vec.end(), v);
    if (it_v == v_vec.end()) {
        // If v does not exist, add it along with the demand value
        v_vec.push_back(v);
        vertex2vertex_value[u].push_back(demand);
    } else {
        // If v already exists, update the demand value
        int index = std::distance(v_vec.begin(), it_v);
        vertex2vertex_value[u][index] = demand;
    }*/
}

void EfficientFRT::removeDemands(double beta, int _v) {
    std::unordered_map<int, std::vector<int>> idVertex2vertices;
    int center = _v;
    for(auto&[u, demand_map] : this->idVertex2idVertex2demand) {
        auto& vertices = idVertex2vertices[u];
        for(auto& [v, d_uv] : demand_map) {
            double firstDistance = graph.getShortestPathDistance(center, u);
            double secondDistance = graph.getShortestPathDistance(center, v);

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



