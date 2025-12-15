//
// Created by Mert Biyikli on 23.10.25.
//

#ifndef OBLIVIOUSROUTING_CKR_PARTITION_H
#define OBLIVIOUSROUTING_CKR_PARTITION_H


#include "../../datastructures/GraphADJ.h"
#include "../../datastructures/GraphCSR.h"
#include "../../datastructures/IGraph.h"
#include "utils/quotient_graph.h"


struct CKRLevel {
    double R = 0.0;                      // the random radius used at this level
    std::vector<int> owner;              // owner[v] = center that captured v at this level (or -1)
    std::vector<int> pred;               // pred[v] = predecessor of v towards its owner center at this level (-1 for center)
    std::vector<int> centers;            // the centers chosen at this level (subset of V)
    std::vector<int> cluster_of_qid; // NEW
};

class CKRPartition {
    GraphADJ m_graph;
public:
    void init(const GraphADJ &g, bool debug=false);

    std::vector<int> computePartition(const std::vector<int>& _X, const double& delta);
    std::vector<int> computePartition(const std::vector<int>& X, const double& delta, CKRLevel& L);
    std::vector<int> computePartition(const IGraph& g, const std::vector<int>& X, const double& delta, CKRLevel& L);
    void build_cluster_of_qid(
    const MendelScaling::QuotientLevel& Q,
    const std::vector<int>& qid_to_rep, // from step 1
    CKRLevel& L
    ) {
        // map center vertex -> cluster index in [0..centers.size)
        std::unordered_map<int,int> center_to_cluster;
        center_to_cluster.reserve(L.centers.size());
        for (int i = 0; i < (int)L.centers.size(); ++i)
            center_to_cluster[L.centers[i]] = i;

        const int nq = (int)qid_to_rep.size();
        L.cluster_of_qid.assign(nq, -1);

        for (int q = 0; q < nq; ++q) {
            int rep = qid_to_rep[q];           // original vertex representative of quotient node q
            int owner_center = L.owner[rep];   // owner of that rep at this CKR level

            auto it = center_to_cluster.find(owner_center);
            if (it == center_to_cluster.end()) {
                throw std::runtime_error("build_cluster_of_qid: owner center not found in centers[]");
            }
            L.cluster_of_qid[q] = it->second;
        }
    }
    // n = #original vertices, nq = #quotient vertices
    std::vector<int> build_qid_to_rep(const MendelScaling::QuotientLevel& Q, int n, int nq) {
        std::vector<int> qid_to_rep(nq, -1);

        // Pick the first original vertex that maps to each quotient id.
        // This is stable if you iterate v in increasing order.
        for (int v = 0; v < n; ++v) {
            int q = Q.sigma_compact_of_v[v];
            if (q < 0 || q >= nq) continue;
            if (qid_to_rep[q] == -1) qid_to_rep[q] = v;
        }

        // Safety: ensure every quotient node got a representative
        for (int q = 0; q < nq; ++q) {
            if (qid_to_rep[q] == -1) {
                throw std::runtime_error("build_qid_to_rep: quotient node without representative");
            }
        }
        return qid_to_rep;
    }

};

#endif //OBLIVIOUSROUTING_CKR_PARTITION_H