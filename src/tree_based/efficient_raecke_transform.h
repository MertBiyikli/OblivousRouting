//
// Created by Mert Biyikli on 24.11.25.
//

#ifndef OBLIVIOUSROUTING_EFFICIENT_RAECKE_TRANSFORM_H
#define OBLIVIOUSROUTING_EFFICIENT_RAECKE_TRANSFORM_H


#include "../graph.h"
#include <unordered_map>
#include <utility>
#include "../utils/hash.h"


/*
 * The efficient Raecke Transform Base class
 * Instead of storing each demand fraction per arc,
 * we store cumulative fractions to optimize the scaling step.
 * Note the following: Since we compute the linear oblivious routing scheme, we can use
 * the linearity by storing each possible commodity s->t as s->x_fixed and t->fixed as in the electrical flow computaiton
 */
template<typename Tree, typename Graph>
class Efficient_RaeckeTransformBase {
protected:
    int x_fixed; // fixed node for demand representation optimization
    std::vector<std::vector<int>> adj_arc_ids; // per-adjacency list version of arc2demand2cumulativeFraction (stores edge ids)
    std::vector<std::vector<double>> adj_arc_demand_cumfractions; // per-adjacency list version of arc2demand2cumulativeFraction

    virtual void normalizeOldSolutionBasedOnNewLambda(double lambda) {
        for (auto& edge_map : adj_arc_demand_cumfractions) {
            for (auto& frac : edge_map) {
                frac *= (1.0 - lambda);
            }
        }
    }

public:
    Efficient_RaeckeTransformBase() {x_fixed = 0;} // default fixed node
    virtual ~Efficient_RaeckeTransformBase() = default;

    virtual void addTree(
        Tree& tree,   // generic pointer, concrete subclasses cast appropriately
        double lambda,
        Graph& g) = 0;
};




#endif //OBLIVIOUSROUTING_EFFICIENT_RAECKE_TRANSFORM_H