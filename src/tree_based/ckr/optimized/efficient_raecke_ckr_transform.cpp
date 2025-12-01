//
// Created by Mert Biyikli on 24.11.25.
//

#include "efficient_raecke_ckr_transform.h"
#include <optional>


void EfficientTreeRoutingTable::init(const int numEdges) {
    adj_ids.resize(numEdges);
    adj_vals.resize(numEdges);
}

std::vector<double> EfficientTreeRoutingTable::operator[](const int e) {
    assert(e >= 0 && e < adj_vals.size());
    return adj_vals[e];
}

void EfficientTreeRoutingTable::addFraction(const int e, const int s, const int t, const double fraction) {
    assert(e >= 0 && e < adj_vals.size());
    auto& ids  = adj_ids[e];
    auto& vals = adj_vals[e];
    int len = static_cast<int>(ids.size());

    // first run linear scan
    int linear_bound = 8;
    for (int i = 0; i < std::min(len, linear_bound); ++i) {
        if (ids[i] == std::make_pair(s, t)) {
            vals[i] += fraction;
            return;
        }
    }

    // binary search to find s, t in ids
    size_t lo = 0, hi = len;
    while (lo < hi) {
        const size_t mid = (lo + hi) >> 1;
        const auto& mid_val = ids[mid];
        if (mid_val < std::make_pair(s, t))
            lo = mid + 1;
        else
            hi = mid;
    }

    if (lo < len && ids[lo] == std::make_pair(s, t)) {
        vals[lo] += fraction;
    } else {
        // Usually append to the end and sort the ids w.r.t. to the terminals
        if (lo == len) {
            ids.emplace_back(s, t);
            vals.push_back(fraction);
        } else {
            ids.insert(ids.begin() + static_cast<long>(lo), {s, t});
            vals.insert(vals.begin() + static_cast<long>(lo), fraction);

        }
    }
}

double EfficientTreeRoutingTable::getFraction(int e , int s, int t) {
    auto& ids  = adj_ids[e];
    auto& vals = adj_vals[e];
    int len = static_cast<int>(ids.size());


    // first run linear scan
    int linear_bound = 8;
    for (int i = 0; i < std::min(len, linear_bound); ++i) {
        if (ids[i] == std::make_pair(s, t)) {
            return vals[i];
        }
    }

    int lo = 0, hi = len;
    while (lo < hi) {
        const size_t mid = (lo + hi) >> 1;
        const auto& mid_val = ids[mid];
        if (mid_val < std::make_pair(s, t))
            lo = mid + 1;
        else
            hi = mid;
    }
    if (lo < len && ids[lo] == std::make_pair(s, t)) {
        return vals[lo];
    } else {
        return 0.0;
    }
}



/*
 * EfficientCKRTransform methods
 */

void EfficientRaeckeCKRTransform::init(IGraph& graph, std::vector<OracleTreeIteration<std::shared_ptr<TreeNode>, std::vector<double> >>& _iters) {
        auto& g_csr = dynamic_cast<Graph_csr&>(graph);
        setGraph(g_csr);
        setIterations(_iters);
    }


void EfficientRaeckeCKRTransform::transform() {
    for (auto& iteration : *iterations) {
        addTree(iteration);
    }
}


void EfficientRaeckeCKRTransform::setGraph(Graph_csr& graph) {
    g = &graph;
    table.init(g->getNumEdges());
}

void EfficientRaeckeCKRTransform::setIterations(std::vector<OracleTreeIteration<std::shared_ptr<TreeNode>, std::vector<double> >>& iterations) {
        this->iterations = &iterations;
}

EfficientTreeRoutingTable& EfficientRaeckeCKRTransform::addTree(OracleTreeIteration<std::shared_ptr<TreeNode>, std::vector<double> >& iteration) {
    distributeDemands(iteration.getTree(), iteration.getLambda(), iteration.getDistance());
    normalizeOldSolutionBasedOnNewLambda(iteration.getLambda());
    removeCycles();
    return table;
}


const EfficientTreeRoutingTable& EfficientRaeckeCKRTransform::getRoutingTable() const {
    return table;
}

void EfficientRaeckeCKRTransform::distributeDemands(const std::shared_ptr<TreeNode> &node, double lambda, const std::vector<double>& distance) {
    if (!node) return;

    for (const auto& child : node->children) {
        std::set<int> A = collectSubtreeVertices(child);
        std::set<int> B;
        for (int v : g->getVertices()) {
            if (!A.count(v)) B.insert(v);
        }

        for (int src : A) {
            for (int dst : B) {
                if (src == dst) continue;
                auto path = g->getShortestPath(src, dst, distance);

                for (size_t i = 0; i + 1 < path.size(); ++i) {
                    std::pair<int, int> arc{path[i], path[i+1]};

                    int edgeId = g->getEdgeId(arc.first, arc.second);
                    int antiEdgeId = g->getEdgeId(arc.second, arc.first);
                    table.addFraction(edgeId, src, dst, lambda);
                    table.addFraction(antiEdgeId, dst, src, lambda);
                }
            }
        }

        distributeDemands(child, lambda, distance);
    }
}

std::set<int> EfficientRaeckeCKRTransform::collectSubtreeVertices(const std::shared_ptr<TreeNode>& node) {
    std::set<int> result(node->members.begin(), node->members.end());
    for (const auto& child : node->children) {
        auto sub = collectSubtreeVertices(child);
        result.insert(sub.begin(), sub.end());
    }
    return result;
}

void EfficientRaeckeCKRTransform::normalizeOldSolutionBasedOnNewLambda(double lambda) {
    for (auto& vals : table.adj_vals) {
        for (auto& frac : vals) {
            frac *= (1.0 - lambda);
        }
    }
}


void EfficientRaeckeCKRTransform::removeCycles() {
    std::set<std::pair<int, int>> all;
    for (int e = 0; e<table.adj_ids.size(); e++) {
        const auto& ids  = table.adj_ids[e];
        for (const auto& id : ids) {
            all.insert(id);
        }
    }

    for (auto d : all) {
        while (true) {
            auto cycle = findCycle(d);
            if (cycle.empty()) break;

            double minF = std::numeric_limits<double>::infinity();
            for (auto e : cycle) {
                int e_id = g->getEdgeId(e.first, e.second);
                double frac = table.getFraction(e_id, d.first, d.second);
                minF = std::min(minF,frac);
            }

            for (auto e : cycle) {
                int e_id = g->getEdgeId(e.first, e.second);
                auto& dmap = table.adj_vals[e_id];
                // find index of (d.first, d.second) in table.adj_ids[e_id]
                const auto& ids = table.adj_ids[e_id];
                size_t len = ids.size();
                size_t lo = 0, hi = len;
                while (lo < hi) {
                    const size_t mid = (lo + hi) >> 1;
                    const auto& mid_val = ids[mid];
                    if (mid_val < std::make_pair(d.first, d.second))
                        lo = mid + 1;
                    else
                        hi = mid;
                }

                // if found, decrease the fraction
                if (lo < len && ids[lo] == std::make_pair(d.first, d.second)) {
                    dmap[lo] -= minF;
                    if (dmap[lo] <= 0) {
                        // remove the entry
                        dmap.erase(dmap.begin() + static_cast<long>(lo));
                        table.adj_ids[e_id].erase(table.adj_ids[e_id].begin() + static_cast<long>(lo));
                    }
                }

                // also erase for the anti edge
                int anti_e_id = g->getEdgeId(e.second, e.first);
                auto& rmap = table.adj_vals[anti_e_id];
                const auto& rids = table.adj_ids[anti_e_id];
                size_t rlen = rids.size();
                size_t rlo = 0, rhi = rlen;
                while (rlo < rhi) {
                    const size_t mid = (rlo + rhi) >> 1;
                    const auto& mid_val = rids[mid];
                    if (mid_val < std::make_pair(d.second, d.first))
                        rlo = mid + 1;
                    else
                        rhi = mid;
                }

                // if found, decrease the fraction
                if (rlo < rlen && rids[rlo] == std::make_pair(d.second, d.first)) {
                    rmap[rlo] -= minF;
                    if (rmap[rlo] <= 0) {
                        // remove the entry
                        rmap.erase(rmap.begin() + static_cast<long>(rlo));
                        table.adj_ids[anti_e_id].erase(table.adj_ids[anti_e_id].begin() + static_cast<long>(rlo));
                    }
                }
            }

        }
    }
}

std::optional<std::vector<std::pair<int, int>>> EfficientRaeckeCKRTransform::findCycleRec(
        int u,
        std::set<int>& analyzed,
        std::vector<int>& stack,
        const std::pair<int, int>& d) {
    if (std::find(stack.begin(), stack.end(), u) != stack.end())
        return std::vector<std::pair<int, int>>{};
    if (analyzed.count(u))
        return std::nullopt;

    analyzed.insert(u);
    stack.push_back(u);

    for (int w : g->neighbors(u)) {
        int e_id = g->getEdgeId(u, w);
        double frac = table.getFraction(e_id, d.first, d.second);
        if (frac <= 0) continue;

        if (auto child = findCycleRec(w, analyzed, stack, d)) {
            auto cycle = *child;
            cycle.insert(cycle.begin(), {u, w});
            return cycle;
        }
    }

    stack.pop_back();
    return std::nullopt;
}

std::vector<std::pair<int, int>> EfficientRaeckeCKRTransform::findCycle(const std::pair<int, int>& d) {
    std::set<int> analyzed;
    for (int v : g->getVertices()) {
        if (analyzed.count(v)) continue;
        std::vector<int> stack;
        if (auto maybe = findCycleRec(v, analyzed, stack, d))
            return *maybe;
    }
    return {};
}