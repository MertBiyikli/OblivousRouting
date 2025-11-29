//
// Created by Mert Biyikli on 27.11.25.
//

#ifndef OBLIVIOUSROUTING_RAECKE_ORACLE_ITERATION_H
#define OBLIVIOUSROUTING_RAECKE_ORACLE_ITERATION_H
#include <memory>


#include <memory>
#include <cassert>

template<typename TreeType, typename DistanceType>
class OracleTreeIteration {
private:
    double lambda = 0.0;
    TreeType tree;          // owns the tree of this level
    DistanceType distance;

public:
    OracleTreeIteration() = default;

    OracleTreeIteration(double lambda,
                        TreeType tree,
                        const DistanceType& distance)
        : lambda(lambda),
          tree((tree)),
          distance(distance)
    {}

    // ---- lambda ----
    void setLambda(double l) noexcept { lambda = l; }
    double getLambda() const noexcept { return lambda; }

    // ---- tree ----
    void setTree(TreeType t) noexcept {
        assert(t != nullptr);
        tree = t;
    }

    const TreeType& getTree() const {
        assert(tree != nullptr);
        return tree;
    }

    TreeType getTree() {
        assert(tree != nullptr);
        return tree;
    }


    // ---- distance ----
    void setDistance(const DistanceType& d) noexcept {
        distance = d;
    }

    const DistanceType& getDistance() const {
        return distance;
    }

};

#endif //OBLIVIOUSROUTING_RAECKE_ORACLE_ITERATION_H