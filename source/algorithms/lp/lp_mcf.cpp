//
// Created by Mert Biyikli on 25.03.26.
//

#include "algorithms/lp/lp_mcf.h"
#include "utils/demands.h"
#include "utils/my_math.h"

Result<void> CMMF_Solver::computeBasisFlows(AllPairRoutingTable& table) {
    // as of now, each LP final class has to call initialize outside the LP base class
    this->n = graph.getNumNodes();
    auto init = this->initSolver();
    if (!init) {
        return getError(init);
    }

    auto var = CreateVariables();
    if (!var) {
        return getError(var);
    }

    CreateConstraints();
    SetObjective();
    // === Solve the LP ===
    status = solver->Solve();
    if (status != MPSolver::OPTIMAL) {
        return makeErrorMessage(ErrorCode::SolverFailed, "Solving the CMF LP returned non-optimal solution.");
    } else {
        storeFlow(table);
    }
    return {};
}


void CMMF_Solver::AddDemandMap(const demands &d_map) {
     for (int i = 0; i < static_cast<int>(d_map.size()); ++i) {
         const auto& d = d_map.getDemandPair(i);
         const double& value = d_map.getDemandValue(i);
         this->addDemand(d.first, d.second, value);
     }
}




Result<void> CMMF_Solver::CreateVariables() {

    alpha = solver->MakeNumVar(0.0, solver->infinity(), "alpha");
    if (!alpha) {
        return makeErrorMessage(ErrorCode::SolverFailed, "Creating the bounding variable for the MCF LP failed.");
    }

    for (int s = 0; s<n; s++) {
        for (int t = 0; t < n; ++t) {
            if (s == t) continue;
            auto &edge2var = map_commodities2edge[{s,t}];
            edge2var.clear();
            for (int e = 0; e<graph.getNumDirectedEdges(); e++) {
                edge2var[e] = solver->MakeNumVar(
                    0.0, solver->infinity(),
                    "f_" + std::to_string(e) + "_" + std::to_string(t)
                );
            }
        }
    }
    return {};
}

void CMMF_Solver::CreateConstraints() {
    // 1) Flow conservation: for each commodity t and node v != t
    //    sum_out f(v->w, t) - sum_in f(w->v, t) = demand(v, t)
    for (int s = 0; s < n; ++s) {
        for (int t = 0; t < n; ++t) {
            if (s == t) continue;

            double rhs = 0.0;
            double it = this->getDemandValue(s, t);
            if (it > EPS  ) rhs = it;

            MPConstraint* const c = solver->MakeRowConstraint(rhs, rhs);

            for (int e = 0; e<graph.getNumDirectedEdges(); e++) {
                const auto &edge = graph.getEdgeEndpoints(e);
                if (edge.first == s)  c->SetCoefficient(map_commodities2edge[{s, t}][e], +1.0); // out
                if (edge.second == s) c->SetCoefficient(map_commodities2edge[{s, t}][e], -1.0); // in
            }
        }
    }

    // 2) Congestion constraints for undirected edges {u,v} handled once
    //    sum_t f(u->v,t) + f(v->u,t) - alpha * cap(u,v) <= 0
    for (int e = 0; e<graph.getNumDirectedEdges(); e++) {
        int u = graph.getEdgeEndpoints(e).first;
        int v = graph.getEdgeEndpoints(e).second;
        if (u > v) continue; // ensure each undirected pair processed once

        // find reverse arc if present
        int rev_id = graph.getAntiEdge(e);

        MPConstraint* const c = solver->MakeRowConstraint(-solver->infinity(), 0.0);

        // -cap(u,v) * alpha
        const double cap = graph.getEdgeCapacity(u, v);
        c->SetCoefficient(alpha, -cap);

        // + sum of flows on both directions (if reverse exists)
        for (int s = 0; s<n; s++) {
            for (int t = 0; t < n; ++t) {
                if (s == t) continue;
                c->SetCoefficient(map_commodities2edge[{s, t}][e], +1.0);
                if (rev_id != -1) {
                    c->SetCoefficient(map_commodities2edge[{s, t}][rev_id], +1.0);
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




void CMMF_Solver::PrintSolution() {
    // 1) Print the congestion factor
    std::cout << "α = " << alpha->solution_value() << "\n\n";

    // 2) For each commodity (i.e. each destination “t”) print all non-zero arc flows
    for (int s = 0; s < n; s++) {
        for (int t = 0; t < n; ++t) {
            if ( s == t) continue;
            std::cout << "Flows for commodity dest=" << s << " -> " << t << ":\n";
            auto &edge2var = map_commodities2edge[{s, t}];
            bool any = false;
            for (auto &kv : edge2var) {
                int arcId    = kv.first;
                MPVariable* v= kv.second;
                double f     = v->solution_value();
                if (f > 1e-8) {
                    any = true;
                    const auto &e = graph.getEdgeEndpoints(arcId);
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
    std::unordered_map<std::pair<int,int>, double, PairHash> und_load;
    for (int s = 0; s < n; s++) {
        for (int t = 0; t < n; ++t) {
            if (s == t) continue;
            for (auto &kv : map_commodities2edge[{s,t}]) {
                int arcId = kv.first;
                double f  = kv.second->solution_value();
                if (f <= SOFT_EPS) continue;
                const auto &e = graph.getEdgeEndpoints(arcId);
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

void CMMF_Solver::storeFlow(AllPairRoutingTable& table) {
    // store the flow for each commodity in f_st_e
    for (int s = 0; s < n; s++) {
        for (int t = 0; t < n; ++t) {
            if (s == t) continue;
            auto &edge2var = map_commodities2edge[{s, t}];
            for (auto &kv : edge2var) {
                int arcId    = kv.first;
                MPVariable* v= kv.second;
                double f     = v->solution_value();
                if (std::abs(f) > SOFT_EPS) {
                    table.addFlow(arcId, s, t, f);
                }
            }
        }
    }
}

Result<double> CMMF_Solver::getCongestionForPassedDemandMap() const {
    if ( solver && alpha
        && status == MPSolver::OPTIMAL) {
        return alpha->solution_value();
    }else {
        return makeErrorMessage(ErrorCode::SolverFailed, "CMMF_Solver: error solving the minimum congestion.");
    }
}
