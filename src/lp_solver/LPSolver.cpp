//
// Created by Mert Biyikli on 11.05.25.
//

#include "LPSolver.h"

using namespace operations_research;


void LPSolver::solve(const Graph& graph) {
    this->Run(graph);
    this->PrintObliviousRoutingSolution(graph.GetDiGraph());
}




void LPSolver::CreateVariables(const DiGraph &graph) {
    if (!solver) solver.reset(MPSolver::CreateSolver("GLOP"));
    if (!solver) throw std::runtime_error("GLOP solver unavailable.");


    // CREATE Alpha variable
    alpha = solver->MakeNumVar(0.0, solver->infinity(), "alpha");

    for(int s = 0; s<graph.numNodes(); s++) {
        for(int t = 0; t<graph.numNodes(); t++) {
            if(s == t ) continue;
            demands.push_back({s, t});
        }
        sourceVertices.push_back(s);
    }


    // p variables for each edge e and i, j
    for( const auto& [id_arc, arc] : graph.GetArcs() ) {
        if (id_arc % 2 == 1) continue;

        for(Demand d : demands) {
            p_e_ij[{id_arc, d.source, d.target}] = solver->MakeNumVar(0.0, solver->infinity(), "p_" + std::to_string(id_arc) + "_" + std::to_string(d.source) + "_" + std::to_string(d.target));
        }


        for(int v = 0; v<graph.numNodes(); v++) {
            p_e_ij[{id_arc, v, v}] = solver->MakeNumVar(0.0, 0.0, "p_" + std::to_string(id_arc) + "_" + std::to_string(v) + "_" + std::to_string(v));
        }
    }





    // π_e_f variables
    for ( const auto& [id_e, e] : graph.GetArcs() ) {
        if( id_e % 2 == 1) continue;

        for ( const auto& [id_f, f] : graph.GetArcs() ) {
            if( id_f % 2 == 1) continue;

            π_e_f[{id_e, id_f}] = solver->MakeNumVar(0.0, solver->infinity(), "w_" + std::to_string(id_e) + "_" + std::to_string(id_f));

        }
    }


    // f_e_st variables
    for(const auto& [id, e] : graph.GetArcs()) {
        for(const auto& d : demands) {
            f_e_st[{id, {d.source, d.target}}] = solver->MakeNumVar(0.0, solver->infinity(), "x_" + std::to_string(id) + "_" + std::to_string(d.source) + "_" + std::to_string(d.target));
            // also t -> s
            // f_e_st[{id, {d.target, d.source}}] = solver->MakeNumVar(0.0, solver->infinity(), "x_" + std::to_string(id)+ "_" + std::to_string(d.target) + "_" + std::to_string(d.source));
        }
    }

}

void LPSolver::CreateConstraints(const DiGraph &graph) {
    // forall links l: sum_{m \in E} of cap(m)*π(l, m) <= alpha
    for(const auto& [id_e, e] : graph.GetArcs()) {
        if ( id_e % 2 == 1) continue;

        MPConstraint* constraint = solver->MakeRowConstraint(-solver->infinity(), 0.0);
        //std::cout << "Dual optimization constraint for arc " << id_e << " and arc " << id_f << ":\n";

        for(const auto& [id_f, f] : graph.GetArcs()) {
            if ( id_f % 2 == 1) continue;

            auto it = π_e_f.find({id_e, id_f});
            if (it != π_e_f.end() && it->second != nullptr) {
                constraint->SetCoefficient(it->second, f.GetCapacity());
                //std::cout << "  1 * π_" << id_e << "^" << id_f << " * " << f.GetCapacity() << "\n";
            }

        }

        constraint->SetCoefficient(alpha, -1);
        //std::cout << "  -1 * alpha * " << "\n";

        // constraint->SetBounds(0, solver->infinity());
    }


    // \forall links l:  f_ij(l)-p_l(i,j)*cap(l) <= 0
    for(auto const& [id, e] : graph.GetArcs()) {
        if ( id % 2 == 1) continue;

        // create a constraint for each arc
        MPConstraint* constraint = solver->MakeRowConstraint(-solver->infinity(), 0.0);
        //std::cout << "Flow constraint for arc " << id << " (" << e.source << " → " << e.target << "):\n";

        for (const auto& d : demands) {
            int s = d.source, t = d.target;

            constraint->SetCoefficient(f_e_st[{id, {s, t}}], 1);
            constraint->SetCoefficient(p_e_ij[{id, s, t}], -e.GetCapacity());


        }

    }

// \forall links l, i, edges e = (j, k) : π(l, link-of-edge(e))+p_l(i,j) - p_l(i,k) >= 0
    for(auto const& [id, e] : graph.GetArcs()) {
        if ( id % 2 == 1) continue;

        for(int i = 0; i<graph.numNodes(); i++) {

            for(const auto&[id_f, f] : graph.GetArcs() ) {
                int j = f.GetSource(), k = f.GetTarget();

                int undirected_link_of_f = (id_f % 2 == 0 ? id_f : id_f - 1);

                MPConstraint* constraint = solver->MakeRowConstraint(0.0, solver->infinity());

                constraint->SetCoefficient(π_e_f[{id, undirected_link_of_f}], 1);

                constraint->SetCoefficient(p_e_ij[{id, i, j}], 1);

                constraint->SetCoefficient(p_e_ij[{id, i, k}], -1);
            }
        }
    }

    // ensure f_ij(e) is a routing

    // TODO: merge the function calls for incoming and outgoing edges together
    // \forall i , j!=i: \sum_{e \in OUT(i)} f_ij(e) - \sum_{e \in IN(i)} f_ij(e) = 1
    for(int i = 0; i<graph.numNodes(); i++) {
        for(int j  = 0; j<graph.numNodes(); j++) {
            if (i == j) continue;

            MPConstraint* constraint = solver->MakeRowConstraint(1.0, 1.0);
            //std::cout << "Routing constraint for demand (" << i << " → " << j << "):\n";

            // outgoing arcs
            for (const auto& [id, e] : graph.GetArcs()) {
                if (e.GetSource() == i) {
                    auto it = f_e_st.find({id, {i, j}});
                    if( it != f_e_st.end() && it->second != nullptr) {
                        constraint->SetCoefficient(it->second, 1);
                    } else {
                        std::cerr << "Warning: Variable for arc " << id << " and demand (" << i << ", " << j << ") not found.\n";
                    }

                    //std::cout << "  +1 * f_" << id << "_(" << i << "," << j << ")\n";
                }
            }

            // incoming arcs
            for (const auto& [id, e] : graph.GetArcs()) {
                if (e.GetTarget() == i) {
                    auto it = f_e_st.find({id, {i, j}});
                    if( it != f_e_st.end() && it->second != nullptr) {
                        constraint->SetCoefficient(it->second, -1);
                    } else {
                        std::cerr << "Warning: Variable for arc " << id << " and demand (" << i << ", " << j << ") not found.\n";
                    }

                    //std::cout << "  -1 * f_" << id << "_(" << i << "," << j << ")\n";
                }
            }
        }
    }

    // \forall k, i != k, j != k: \sum_{e \in OUT(k)} f_ij(e) - \sum_{e \in IN(k)} f_ij(e) = 0
    for(int k = 0; k<graph.numNodes(); k++) {
        for(int i = 0; i<graph.numNodes(); i++) {
            if (i == k) continue;
            for(int j = 0; j<graph.numNodes(); j++) {
                if (j == k || j == i) continue;

                MPConstraint* constraint = solver->MakeRowConstraint(0.0, 0.0);
                //std::cout << "Flow conservation constraint at node " << k << " for demand (" << i << " → " << j << "):\n";

                // outgoing arcs
                for (const auto& [id, e] : graph.GetArcs()) {
                    if (e.GetSource() == k) {
                        auto it = f_e_st.find({id, {i, j}});
                        if (it != f_e_st.end() && it->second != nullptr) {
                            constraint->SetCoefficient(it->second, 1);
                        } else {
                            std::cerr << "Warning: Variable for arc " << id << " and demand (" << i << ", " << j << ") not found.\n";
                        }

                        //std::cout << "  +1 * f_" << id << "_(" << i << "," << j << ")\n";
                    }
                }

                // incoming arcs
                for (const auto& [id, e] : graph.GetArcs()) {
                    if (e.GetTarget() == k) {
                        auto it = f_e_st.find({id, {i, j}});
                        if (it != f_e_st.end() && it->second != nullptr) {
                            constraint->SetCoefficient(it->second, -1);
                        } else {
                            std::cerr << "Warning: Variable for arc " << id << " and demand (" << i << ", " << j << ") not found.\n";
                        }

                        //std::cout << "  -1 * f_" << id << "_(" << i << "," << j << ")\n";
                    }
                }
            }
        }
    }

}

void LPSolver::GetRoutingTable(const DiGraph& graph) {
    std::cout << "\n=== Oblivious Routing Table ===\n";

    for (const auto& [key, var] : f_e_st) {

        if (!var) {
            std::cerr << "Warning: variable pointer is null for key.\n";
            continue;
        }
        double value = var->solution_value();
        const auto& arc = std::get<0>(key);
        const std::pair<int, int>& demand = std::get<1>(key);



        if (value > 1e-6) { // skip near-zero
            std::cout << "Arc " << arc << " (" << graph.GetArcById(arc).GetSource() << " → " << graph.GetArcById(arc).GetTarget() << ") "
                      << "carries " << value << " units of demand "
                      << "from " << demand.first << " to " << demand.second << "\n";
        }
    }
}

void LPSolver::PrintObliviousRoutingSolution(const DiGraph &graph) {

    // === Print the solution ===
    std::unordered_map<int, double> total_flow_per_arc;

    for (const auto& [key, var] : f_e_st) {
        if (!var) continue;
        double val = var->solution_value();
        int arc_id = std::get<0>(key);
        total_flow_per_arc[arc_id] += val;
    }

    for (const auto& [id, arc] : graph.GetArcs()) {
        if (id % 2 == 1) continue;

        double total_flow = total_flow_per_arc[id];
        if(max_cong < total_flow/arc.GetCapacity()) {
            max_cong = total_flow/arc.GetCapacity();
        }

    }
    std::cout << "Max congestion across all undirected arcs: " << max_cong << std::endl;

}
