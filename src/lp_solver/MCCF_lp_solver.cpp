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


void CMMF_Solver::CreateVariables(const RaeckeGraph &graph) {


    // CREATE Alpha variable
    alpha = solver->MakeNumVar(0.0, solver->infinity(), "alpha");

    for(int v = 0; v<n; v++) {
        this->map_vertex2edge[v] = std::unordered_map<int, MPVariable*>();
        for(int edge_id = 0;  edge_id<edges.size(); edge_id++) {
            if(edges[edge_id].first == v) continue;

            this->map_vertex2edge[v][edge_id] =solver->MakeNumVar(0.0, solver->infinity(),
                                                                  "f_" + std::to_string(edge_id) + "_" + std::to_string(v));
        }
    }
}

void CMMF_Solver::CreateConstraints(const RaeckeGraph &graph) {

    for (int dest = 0; dest < n; ++dest) {
        for (int v = 0; v < n; ++v) {
            if (v == dest) continue;

            double rhs = 0.0;
            auto demand_it = m_demands.find({v, dest});
            if(demand_it != m_demands.end()) rhs = demand_it->second;

            MPConstraint* constraint = solver->MakeRowConstraint(rhs, solver->infinity());


            for(int edge_id = 0; edge_id<edges.size(); edge_id++) {

                if(edges[edge_id].first == v) {
                    constraint->SetCoefficient(map_vertex2edge[dest][edge_id], 1.0);

                    // get reverse edge
                    auto rev_it = std::find(edges.begin(), edges.end(), std::make_pair(edges[edge_id].second, edges[edge_id].first));
                    if(rev_it != edges.end()
                    && ((*rev_it).first != dest)
                    ) {
                        int rev_edge_id = std::distance(edges.begin(),rev_it);
                        constraint->SetCoefficient( map_vertex2edge[dest][rev_edge_id], -1.0);
                    }
                }
            }

            // add constraint greater or equal

        }
    }


    for(int edge_id = 0; edge_id<edges.size(); edge_id++) {
        int u = edges[edge_id].first, v = edges[edge_id].second;
        if(u > v) continue;

        MPConstraint* constraint = solver->MakeRowConstraint(0.0, solver->infinity());

        constraint->SetCoefficient(alpha,  graph.getEdgeCapacity(u, v));

        for(int dest = 0; dest < n; dest++) {
            if(u == dest) continue;
            constraint->SetCoefficient(map_vertex2edge[dest][edge_id], -1.0);
            if(v == dest) continue;

            auto rev_it = std::find(edges.begin(), edges.end(), std::make_pair(v,v));
            int rev_edge_id = -1;
            if(rev_it != edges.end()) {
                rev_edge_id = std::distance(edges.begin(), rev_it);
                constraint->SetCoefficient(map_vertex2edge[dest][rev_edge_id], -1.0);
            }
        }
    }

}

void CMMF_Solver::SetObjective() {
    // === OBJECTIVE ===
    solver->MutableObjective()->SetCoefficient(alpha, 1.0);
    solver->MutableObjective()->SetMinimization();
}

void CMMF_Solver::PrintSolution(const RaeckeGraph &graph) {
    // 1) Print the congestion factor
    std::cout << "α = " << alpha->solution_value() << "\n\n";

    // 2) For each commodity (i.e. each destination “t”) print all non-zero arc flows
    for (int t = 0; t < graph.getNumNodes(); ++t) {
        std::cout << "Flows for commodity dest=" << t << ":\n";
        auto &edge2var = map_vertex2edge[t];
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

    // 3) (Optional) If you also want to see total load per undirected edge:
    //    sum over both directions
    std::cout << "Total load per undirected edge:\n";
    std::map<std::pair<int,int>, double> und_load;
    for (int t = 0; t < graph.getNumNodes(); ++t) {
        for (auto &kv : map_vertex2edge[t]) {
            int arcId = kv.first;
            double f  = kv.second->solution_value();
            if (f <= 1e-8) continue;
            const auto &e = edges[arcId];
            // key = sorted pair of endpoints
            auto key = std::minmax(e.first, e.second);
            und_load[key] += f;
        }
    }
    for (auto &kv : und_load) {
        auto [u,v] = kv.first;
        std::cout
                << "  {" << u << "," << v << "}: "
                << kv.second << "\n";
    }
}
