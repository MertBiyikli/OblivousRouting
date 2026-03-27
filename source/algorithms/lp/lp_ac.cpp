//
// Created by Mert Biyikli on 25.03.26.
//

#include "../../../include/algorithms/lp/lp_ac.h"
#include "../../../include/utils/hash.h"
#include "../../../include/utils/my_math.h"


using namespace operations_research;


void LPSolver::CreateVariables() {
    if (!solver) solver.reset(MPSolver::CreateSolver("GLOP"));
    if (!solver) throw std::runtime_error("GLOP solver unavailable.");


    // CREATE Alpha variable
    alpha = solver->MakeNumVar(0.0, solver->infinity(), "alpha");

    for(int s = 0; s<n; s++) {
        for(int t = 0; t<n; t++) {
            if(s == t ) continue;
            m_demands.push_back({s, t});
        }
        sourceVertices.push_back(s);
    }


    // p variables for each edge e and i, j
    for(int e = 0; e < graph.getNumDirectedEdges(); e++) {
        const auto [u, v] = graph.getEdgeEndpoints(e);
        if (u > v) continue;
        for(auto d : m_demands) {
            p_e_ij[{e, d.first, d.second}] = solver->MakeNumVar(0.0, solver->infinity(), "p_" + std::to_string(e) + "_" + std::to_string(d.first) + "_" + std::to_string(d.second));
        }
        for(int w = 0; w<n; w++) {
            p_e_ij[{e, w, w}] = solver->MakeNumVar(0.0, 0.0, "p_" + std::to_string(e) + "_" + std::to_string(w) + "_" + std::to_string(w));
        }
    }





    // π_e_f variables
    for(int e = 0; e < graph.getNumDirectedEdges(); e++) {
        const auto& [u, v] = graph.getEdgeEndpoints(e);
        if (u > v) continue;

        for ( int f = 0; f < graph.getNumDirectedEdges();f++) {
            if(graph.getEdgeEndpoints(f).first > graph.getEdgeEndpoints(f).second) continue;

            π_e_f[{e, f}] = solver->MakeNumVar(0.0, solver->infinity(), "w_" + std::to_string(e) + "_" + std::to_string(f));

        }
    }


    // f_e_st variables for all arcs
    for(int e = 0; e < graph.getNumDirectedEdges(); e++) {
        for(const auto& d : m_demands) {
            m_var_f_e_[{e, {d.first, d.second}}] = solver->MakeNumVar(0.0, solver->infinity(), "x_" + std::to_string(e) + "_" + std::to_string(d.first) + "_" + std::to_string(d.second));
            // also t -> s
            // f_e_st[{id, {d.target, d.source}}] = solver->MakeNumVar(0.0, solver->infinity(), "x_" + std::to_string(id)+ "_" + std::to_string(d.target) + "_" + std::to_string(d.source));
        }
    }

}

void LPSolver::CreateConstraints() {
    // forall links l: sum_{m \in E} of cap(m)*π(l, m) <= alpha
    for(int e = 0; e < graph.getNumDirectedEdges(); e++) {
        if ( graph.getEdgeEndpoints(e).first > graph.getEdgeEndpoints(e).second) continue;

        MPConstraint* constraint = solver->MakeRowConstraint(-solver->infinity(), 0.0);
        //std::cout << "Dual optimization constraint for arc " << id_e << " and arc " << id_f << ":\n";

        for(int f = 0; f < graph.getNumDirectedEdges(); f++) {
            if ( graph.getEdgeEndpoints(f).first > graph.getEdgeEndpoints(f).second ) continue;

            auto it = π_e_f.find({e, f});
            if (it != π_e_f.end()
                && it->second) {
                constraint->SetCoefficient(it->second, graph.getEdgeCapacity(e));
            }
        }

        constraint->SetCoefficient(alpha, -1);
    }


    // \forall links l:  f_ij(l)-p_l(i,j)*cap(l) <= 0
    for(int e = 0; e < graph.getNumDirectedEdges(); e++) {
        if ( graph.getEdgeEndpoints(e).first > graph.getEdgeEndpoints(e).second) continue;

        // create a constraint for each arc
        MPConstraint* constraint = solver->MakeRowConstraint(-solver->infinity(), 0.0);

        for (const auto& d : m_demands) {
            int s = d.first, t = d.second;

            constraint->SetCoefficient(m_var_f_e_[{e, {s, t}}], 1);
            constraint->SetCoefficient(p_e_ij[{e, s, t}], -graph.getEdgeCapacity(e));


        }

    }

// \forall links l, i, edges e = (j, k) : π(l, link-of-edge(e))+p_l(i,j) - p_l(i,k) >= 0
    for(int e = 0; e < graph.getNumDirectedEdges(); e++) {
        if ( graph.getEdgeEndpoints(e).first > graph.getEdgeEndpoints(e).second) continue;

        for(int i = 0; i<n; i++) {

            for( int f = 0; f < graph.getNumDirectedEdges(); f++) {
                int j = graph.getEdgeEndpoints(f).first, k = graph.getEdgeEndpoints(f).second;

                int undirected_link_of_f = -1;
                if (j > k) {
                    undirected_link_of_f = graph.getAntiEdge(f);
                }else {
                    undirected_link_of_f = f;
                }


                MPConstraint* constraint = solver->MakeRowConstraint(0.0, solver->infinity());
                constraint->SetCoefficient(π_e_f[{e, undirected_link_of_f}], 1);
                constraint->SetCoefficient(p_e_ij[{e, i, j}], 1);
                constraint->SetCoefficient(p_e_ij[{e, i, k}], -1);
            }
        }
    }

    // \forall i , j!=i: \sum_{e \in OUT(i)} f_ij(e) - \sum_{e \in IN(i)} f_ij(e) = 1
    for(int i = 0; i<n; i++) {
        for(int j  = 0; j<n; j++) {
            if (i == j) continue;

            MPConstraint* constraint = solver->MakeRowConstraint(1.0, 1.0);

            // outgoing arcs
            for ( int e = 0; e < graph.getNumDirectedEdges(); e++ ) {

                if (graph.getEdgeEndpoints(e).first == i) {
                    auto it = m_var_f_e_.find({e, {i, j}});
                    if( it != m_var_f_e_.end() &&
                        it->second) {
                        constraint->SetCoefficient(it->second, 1);
                    } else {
                        std::cerr << "Warning: Variable for arc " << e << " and demand (" << i << ", " << j << ") not found.\n";
                    }
                }
            }

            // incoming arcs
            for (int e = 0; e < graph.getNumDirectedEdges(); e++ ) {
                if (graph.getEdgeEndpoints(e).second == i) {
                    auto it = m_var_f_e_.find({e, {i, j}});
                    if( it != m_var_f_e_.end()
                        && it->second) {
                        constraint->SetCoefficient(it->second, -1);
                    } else {
                        std::cerr << "Warning: Variable for arc " << e << " and demand (" << i << ", " << j << ") not found.\n";
                    }
                }
            }
        }
    }

    // \forall k, i != k, j != k: \sum_{e \in OUT(k)} f_ij(e) - \sum_{e \in IN(k)} f_ij(e) = 0
    for(int k = 0; k<n; k++) {
        for(int i = 0; i<n; i++) {
            if (i == k) continue;
            for(int j = 0; j<n; j++) {
                if (j == k || j == i) continue;

                MPConstraint* constraint = solver->MakeRowConstraint(0.0, 0.0);
                // outgoing arcs
                for (int e = 0; e < graph.getNumDirectedEdges(); e++ ) {
                    if (graph.getEdgeEndpoints(e).first == k) {
                        auto it = m_var_f_e_.find({e, {i, j}});
                        if (it != m_var_f_e_.end() && it->second != nullptr) {
                            constraint->SetCoefficient(it->second, 1);
                        } else {
                            std::cerr << "Warning: Variable for arc " << e << " and demand (" << i << ", " << j << ") not found.\n";
                        }
                    }
                }

                // incoming arcs
                for (int e = 0; e < graph.getNumDirectedEdges(); e++) {
                    if (graph.getEdgeEndpoints(e).second == k) {
                        auto it = m_var_f_e_.find({e, {i, j}});
                        if (it != m_var_f_e_.end() && it->second != nullptr) {
                            constraint->SetCoefficient(it->second, -1);
                        } else {
                            std::cerr << "Warning: Variable for arc " << e << " and demand (" << i << ", " << j << ") not found.\n";
                        }
                    }
                }
            }
        }
    }
}

void LPSolver::SetObjective()
{
    // === Objective: maximize alpha ===
    solver->MutableObjective()->SetCoefficient(alpha, 1);
    solver->MutableObjective()->SetMinimization();
}



void LPSolver::PrintSolution() {

    max_cong = 0;
    // === Print the solution ===
    std::unordered_map<int, double> total_flow_per_arc;

    for (const auto& [key, var] : m_var_f_e_) {
        if (!var) continue;
        double val = var->solution_value();
        int arc_id = std::get<0>(key);
        total_flow_per_arc[arc_id] += val;
    }

    for (int e = 0; e < graph.getNumDirectedEdges(); e++) {
        int u = graph.getEdgeEndpoints(e).first, v = graph.getEdgeEndpoints(e).second;

        if (u > v) continue;

        int rev_edge = graph.getAntiEdge(e);
        double total_flow = total_flow_per_arc[e];

        // to ensure that flow along anti-parallel arcs
        // will be added as absolute flow to its corresponding undirected link
        total_flow += total_flow_per_arc[rev_edge];

        double capacity = graph.getEdgeCapacity(e);

        if(max_cong < total_flow/capacity) {
            max_cong = total_flow/capacity;
        }

    }
    std::cout << "Max congestion across all undirected arcs: " << max_cong << std::endl;

}



void LPSolver::storeFlow(AllPairRoutingTable& table) {

    for (int e = 0; e<graph.getNumDirectedEdges(); e++) {

        // Iterate over all demands for this edge
        for (const auto& d : m_demands) {
            auto it = m_var_f_e_.find({e, {d.first, d.second}});


            if (it != m_var_f_e_.end()
                && it->second) {

                double flow_value = it->second->solution_value();
                if (std::abs(flow_value) > SOFT_EPS) { // print only non-zero flows
                    //int e = graph.getEdgeId(graph.getEdgeEndpoints(e).first, graph.getEdgeEndpoints(e).second);

                    if (flow_value < 0) {
                        // push flow into anti-edge direction
                        int anti_e = graph.getAntiEdge(e);
                        table.addFlow(anti_e, d.first, d.second, flow_value);
                    }else {
                        table.addFlow(e, d.first, d.second, flow_value);
                    }
                }
            }
        }
    }
}
