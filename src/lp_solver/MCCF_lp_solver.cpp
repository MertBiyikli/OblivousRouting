//
// Created by Mert Biyikli on 30.06.25.
//

#include "MCCF_lp_solver.h"
#include "../utils/hash.h"




void CMMF_Solver::AddDemands(const Demand& d, double demand) {
    int u = d.source;
    int v = d.target;
    if (u < 0 || v < 0
        || u >= n || v >= n) {
        throw std::invalid_argument("Vertex IDs out of range.");
    }
    if (u == v) {
        throw std::invalid_argument("Cannot add demand from a vertex to itself.");
    }

    if(m_demands.find(d) != m_demands.end()) {
        throw std::invalid_argument("Demand between these vertices already exists.");
    }


    m_demands[d] = demand;
}


// C++
void CMMF_Solver::CreateVariables(const Graph &graph) {
    // α: maximum congestion
    alpha = solver->MakeNumVar(0.0, solver->infinity(), "alpha");

    // f[eid, t] >= 0 for every undirected edge and destination t
    for (int s = 0; s<n; s++) {
        for (int t = 0; t < n; ++t) {
            if (s == t) continue;
            auto &edge2var = map_vertex2edge[{s,t}];
            edge2var.clear();
            for (int eid = 0; eid < static_cast<int>(edges.size()); ++eid) {
                edge2var[eid] = solver->MakeNumVar(
                    0.0, solver->infinity(),
                    "f_" + std::to_string(eid) + "_" + std::to_string(t)
                );
            }
        }
    }
}

void CMMF_Solver::CreateConstraints(const Graph &graph) {
    // 1) Flow conservation: for each commodity t and node v != t
    //    sum_out f(v->w, t) - sum_in f(w->v, t) = demand(v, t)
    for (int s = 0; s < n; ++s) {
        for (int t = 0; t < n; ++t) {
            if (s == t) continue;

            double rhs = 0.0;
            auto it = m_demands.find({s, t});
            if (it != m_demands.end()) rhs = it->second;

            MPConstraint* c = solver->MakeRowConstraint(rhs, rhs);

            for (int eid = 0; eid < static_cast<int>(edges.size()); ++eid) {
                const auto &e = edges[eid];
                if (e.first == s)  c->SetCoefficient(map_vertex2edge[{s, t}][eid], +1.0); // out
                if (e.second == s) c->SetCoefficient(map_vertex2edge[{s, t}][eid], -1.0); // in
            }
        }
    }

    // 2) Congestion (capacity) constraints for undirected edges {u,v} handled once
    //    sum_t f(u->v,t) + f(v->u,t) - alpha * cap(u,v) <= 0
    for (int eid = 0; eid < static_cast<int>(edges.size()); ++eid) {
        int u = edges[eid].first;
        int v = edges[eid].second;
        if (u > v) continue; // ensure each undirected pair processed once

        // find reverse arc if present
        int rev_id = -1;
        auto rev_it = std::find(edges.begin(), edges.end(), std::make_pair(v, u));
        if (rev_it != edges.end()) rev_id = static_cast<int>(std::distance(edges.begin(), rev_it));

        MPConstraint* c = solver->MakeRowConstraint(-solver->infinity(), 0.0);

        // -cap(u,v) * alpha
        const double cap = graph.getEdgeCapacity(u, v);
        c->SetCoefficient(alpha, -cap);

        // + sum of flows on both directions (if reverse exists)
        for (int s = 0; s<n; s++) {
            for (int t = 0; t < n; ++t) {
                if (s == t) continue;
                c->SetCoefficient(map_vertex2edge[{s, t}][eid], +1.0);
                if (rev_id != -1) {
                    c->SetCoefficient(map_vertex2edge[{s, t}][rev_id], +1.0);
                }
            }
        }
    }
}

void CMMF_Solver::SetObjective() {
    auto* obj = solver->MutableObjective();
    obj->SetCoefficient(alpha, 1.0);
    obj->SetMinimization();
}


void CMMF_Solver::PrintSolution(const Graph &graph) {
    // 1) Print the congestion factor
    std::cout << "α = " << alpha->solution_value() << "\n\n";

    // 2) For each commodity (i.e. each destination “t”) print all non-zero arc flows
    for (int s = 0; s < n; s++) {
        for (int t = 0; t < graph.getNumNodes(); ++t) {
            if ( s == t) continue;
            std::cout << "Flows for commodity dest=" << s << " -> " << t << ":\n";
            auto &edge2var = map_vertex2edge[{s, t}];
            bool any = false;
            for (auto &kv : edge2var) {
                int arcId    = kv.first;
                MPVariable* v= kv.second;
                double f     = v->solution_value();
                if (f > 1e-8) {
                    any = true;
                    const auto &e = edges[arcId];
                    std::cout
                            << "  edge (" << e.first << " / " << e.second << "): "
                            << f << "\n";
                }
            }
            if (!any) {
                std::cout << "  (all zero)\n";
            }
            std::cout << "\n";
        }
    }
    // 3) (Optional) If you also want to see total load per undirected edge:
    //    sum over both directions
    std::cout << "Total load per undirected edge:\n";
    std::map<std::pair<int,int>, double> und_load;
    for (int s = 0; s < n; s++) {
        for (int t = 0; t < graph.getNumNodes(); ++t) {
            if (s == t) continue;
            for (auto &kv : map_vertex2edge[{s,t}]) {
                int arcId = kv.first;
                double f  = kv.second->solution_value();
                if (f <= 1e-8) continue;
                const auto &e = edges[arcId];
                // key = sorted pair of endpoints
                auto key = std::minmax(e.first, e.second);
                und_load[key] += f;
            }
        }
    }
    for (auto &kv : und_load) {
        auto [u,v] = kv.first;
        std::cout
                << "  {" << u << "," << v << "}: "
                << kv.second << "\n";
    }
}

void CMMF_Solver::storeFlow() {

}