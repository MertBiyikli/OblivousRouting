//
// Created by Mert Biyikli on 30.06.25.
//

#include "MCCF_lp_solver.h"
#include "../utils/hash.h"

void CMMF_Solver::solve(const Graph &graph) {
    this->Run(graph);
    this->PrintSolution(graph.GetDiGraph());
}

void CMMF_Solver::CreateVariables(const DiGraph &graph) {
    if (!solver) solver.reset(MPSolver::CreateSolver("GLOP"));
    if (!solver) throw std::runtime_error("GLOP solver unavailable.");

    // CREATE Alpha variable
    alpha = solver->MakeNumVar(0.0, solver->infinity(), "alpha");


    for(int s = 0; s < graph.numNodes(); s++) {
        for (const auto &[id, e]: graph.GetArcs()) {

            if (e.source == s) continue;

            map_vertex2edge[s][id] = solver->MakeNumVar(0.0, solver->infinity(),
                                                        "map_" + std::to_string(s) + "_" + std::to_string(id));
        }
    }

}

void CMMF_Solver::AddDemands(const Demand& d, double demand) {
    int u = d.source;
    int v = d.target;
    if (u < 0 || v < 0) {
        throw std::invalid_argument("Vertex IDs must be non-negative.");
    }
    if (u == v) {
        throw std::invalid_argument("Cannot add demand from a vertex to itself.");
    }

    if(m_demands.find(d) != m_demands.end()) {
        throw std::invalid_argument("Demand between these vertices already exists.");
    }
    m_demands[d] = demand;
}

void CMMF_Solver::CreateConstraints(const DiGraph &graph) {


    // Build quick access to incoming arcs
    // incoming[v] = list of arc‐ids whose head is v
    std::vector<std::vector<int>> incoming(graph.numNodes());
    for (auto const& [id, e] : graph.GetArcs()) {
        incoming[e.target].push_back(id);
    }

    // Define flow conservation constraints
    for(int t = 0; t < graph.numNodes(); t++) {
        for(int v = 0; v < graph.numNodes(); v++) {

            if(v == t) continue; // Skip if source and target are the same

            double demand_value = 0.0;
            auto it = m_demands.find({v, t});
            if(it != m_demands.end()) {
                demand_value = it->second;
            }

            MPConstraint* constraint = solver->MakeRowConstraint(demand_value, solver->infinity());


            // outgoing neighbors
            for(const auto& [id, arc] : graph.GetArcs()) {
                if(arc.source != v) continue; // Only consider outgoing arcs from v

                constraint->SetCoefficient(map_vertex2edge[t][arc.GetId()], 1.0);

                auto rev_arc = (arc.GetReverseArcId());

                if(arc.target != t) {
                    constraint->SetCoefficient(map_vertex2edge[t][rev_arc], -1.0);
                }
            }

        }
    }

    for (auto const& [arcId, e] : graph.GetArcs()) {
        if(arcId % 2 == 1) continue; // Skip odd arc IDs, as per the original code logic

        // LHS: α·C_uv  − ∑_t g_t(u→v)  ≥  0
        MPConstraint* c = solver->MakeRowConstraint(0.0, solver->infinity(),
                                                    "cap_" + std::to_string(arcId));

        c->SetCoefficient(alpha, e.capacity);

        for(int v(0); v<graph.numNodes(); v++) {

            if(e.source == v) continue;

            c->SetCoefficient(map_vertex2edge[v][arcId], -1.0); // Coefficient for the edge variable


            if(e.target == v) continue; // Skip if the vertex is the target

            c->SetCoefficient(map_vertex2edge[v][e.GetReverseArcId()], -1.0); // Coefficient for the reverse edge variable

        }
    }

}

void CMMF_Solver::PrintSolution(const DiGraph &graph) {
    // 1) Print the congestion factor
    std::cout << "α = " << alpha->solution_value() << "\n\n";

    // 2) For each commodity (i.e. each destination “t”) print all non-zero arc flows
    for (int t = 0; t < graph.numNodes(); ++t) {
        std::cout << "Flows for commodity dest=" << t << ":\n";
        auto &edge2var = map_vertex2edge[t];
        bool any = false;
        for (auto &kv : edge2var) {
            int arcId    = kv.first;
            MPVariable* v= kv.second;
            double f     = v->solution_value();
            if (f > 1e-8) {
                any = true;
                const auto &e = graph.GetArcById(arcId);
                std::cout
                        << "  arc " << arcId
                        << " (" << e.source << "→" << e.target << "): "
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
    for (int t = 0; t < graph.numNodes(); ++t) {
        for (auto &kv : map_vertex2edge[t]) {
            int arcId = kv.first;
            double f  = kv.second->solution_value();
            if (f <= 1e-8) continue;
            const auto &e = graph.GetArcById(arcId);
            // key = sorted pair of endpoints
            auto key = std::minmax(e.source, e.target);
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

