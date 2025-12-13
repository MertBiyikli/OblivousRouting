//
// Created by Mert Biyikli on 11.05.25.
//

#include "LPSolver.h"

using namespace operations_research;


void LPSolver::CreateVariables(const IGraph &graph) {
    if (!solver) solver.reset(MPSolver::CreateSolver("GLOP"));
    if (!solver) throw std::runtime_error("GLOP solver unavailable.");


    // CREATE Alpha variable
    alpha = solver->MakeNumVar(0.0, solver->infinity(), "alpha");

    for(int s = 0; s<graph.getNumNodes(); s++) {
        for(int t = 0; t<graph.getNumNodes(); t++) {
            if(s == t ) continue;
            demands.push_back({s, t});
        }
        sourceVertices.push_back(s);
    }


    // p variables for each edge e and i, j
    for(int id = 0; id<edges.size() ; id++ ) {
        if (edges[id].first > edges[id].second) continue;

        for(Demand d : demands) {
            p_e_ij[{id, d.first, d.second}] = solver->MakeNumVar(0.0, solver->infinity(), "p_" + std::to_string(id) + "_" + std::to_string(d.first) + "_" + std::to_string(d.second));
        }


        for(int v = 0; v<graph.getNumNodes(); v++) {
            p_e_ij[{id, v, v}] = solver->MakeNumVar(0.0, 0.0, "p_" + std::to_string(id) + "_" + std::to_string(v) + "_" + std::to_string(v));
        }
    }





    // π_e_f variables
    for ( int id_e = 0; id_e < edges.size(); id_e++) {
        if(edges[ id_e ].first > edges[id_e].second) continue;

        for ( int id_f = 0; id_f < edges.size(); id_f++ ) {
            if( edges[id_f].first > edges[id_f].second) continue;

            π_e_f[{id_e, id_f}] = solver->MakeNumVar(0.0, solver->infinity(), "w_" + std::to_string(id_e) + "_" + std::to_string(id_f));

        }
    }


    // f_e_st variables
    for(int id = 0; id < edges.size(); id++) {
        for(const auto& d : demands) {
            m_var_f_e_[{id, {d.first, d.second}}] = solver->MakeNumVar(0.0, solver->infinity(), "x_" + std::to_string(id) + "_" + std::to_string(d.first) + "_" + std::to_string(d.second));
            // also t -> s
            // f_e_st[{id, {d.target, d.source}}] = solver->MakeNumVar(0.0, solver->infinity(), "x_" + std::to_string(id)+ "_" + std::to_string(d.target) + "_" + std::to_string(d.source));
        }
    }

}

void LPSolver::CreateConstraints(const IGraph &graph) {
    // forall links l: sum_{m \in E} of cap(m)*π(l, m) <= alpha
    for(int id = 0; id<edges.size(); id++) {
        if ( edges[id].first > edges[id].second) continue;

        MPConstraint* constraint = solver->MakeRowConstraint(-solver->infinity(), 0.0);
        //std::cout << "Dual optimization constraint for arc " << id_e << " and arc " << id_f << ":\n";

        for(int id_f = 0; id_f < edges.size(); id_f++) {
            if ( edges[id].first > edges[id].second ) continue;

            auto it = π_e_f.find({id, id_f});
            if (it != π_e_f.end() && it->second != nullptr) {
                constraint->SetCoefficient(it->second, graph.getEdgeCapacity(edges[id].first, edges[id].second));
                //std::cout << "  1 * π_" << id_e << "^" << id_f << " * " << f.GetCapacity() << "\n";
            }

        }

        constraint->SetCoefficient(alpha, -1);
        //std::cout << "  -1 * alpha * " << "\n";

        // constraint->SetBounds(0, solver->infinity());
    }


    // \forall links l:  f_ij(l)-p_l(i,j)*cap(l) <= 0
    for(int id = 0; id< edges.size(); id++) {
        if ( edges[id].first > edges[id].second) continue;

        // create a constraint for each arc
        MPConstraint* constraint = solver->MakeRowConstraint(-solver->infinity(), 0.0);
        //std::cout << "Flow constraint for arc " << id << " (" << e.source << " → " << e.target << "):\n";

        for (const auto& d : demands) {
            int s = d.first, t = d.second;

            constraint->SetCoefficient(m_var_f_e_[{id, {s, t}}], 1);
            constraint->SetCoefficient(p_e_ij[{id, s, t}], -graph.getEdgeCapacity(edges[id].first, edges[id].second));


        }

    }

// \forall links l, i, edges e = (j, k) : π(l, link-of-edge(e))+p_l(i,j) - p_l(i,k) >= 0
    for(int id = 0; id<edges.size(); id++) {
        if ( edges[id].first > edges[id].second) continue;

        for(int i = 0; i<graph.getNumNodes(); i++) {

            for( int id_f = 0; id_f < edges.size(); id_f++ ) {
                int j = edges[id_f].first, k = edges[id_f].second;

                int undirected_link_of_f = -1;

                auto undirected_edge_it = std::find(edges.begin(), edges.end(), std::make_pair(k, j));
                if(undirected_edge_it != edges.end()) {
                    undirected_link_of_f = std::distance(edges.begin(),undirected_edge_it);
                }


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
    for(int i = 0; i<graph.getNumNodes(); i++) {
        for(int j  = 0; j<graph.getNumNodes(); j++) {
            if (i == j) continue;

            MPConstraint* constraint = solver->MakeRowConstraint(1.0, 1.0);
            //std::cout << "Routing constraint for demand (" << i << " → " << j << "):\n";

            // outgoing arcs
            for ( int id = 0; id<edges.size(); id++ ) {

                if (edges[id].first == i) {
                    auto it = m_var_f_e_.find({id, {i, j}});
                    if( it != m_var_f_e_.end() && it->second != nullptr) {
                        constraint->SetCoefficient(it->second, 1);
                    } else {
                        std::cerr << "Warning: Variable for arc " << id << " and demand (" << i << ", " << j << ") not found.\n";
                    }

                    //std::cout << "  +1 * f_" << id << "_(" << i << "," << j << ")\n";
                }
            }

            // incoming arcs
            for (int id = 0; id<edges.size(); id++) {
                if (edges[id].second == i) {
                    auto it = m_var_f_e_.find({id, {i, j}});
                    if( it != m_var_f_e_.end() && it->second != nullptr) {
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
    for(int k = 0; k<graph.getNumNodes(); k++) {
        for(int i = 0; i<graph.getNumNodes(); i++) {
            if (i == k) continue;
            for(int j = 0; j<graph.getNumNodes(); j++) {
                if (j == k || j == i) continue;

                MPConstraint* constraint = solver->MakeRowConstraint(0.0, 0.0);
                //std::cout << "Flow conservation constraint at node " << k << " for demand (" << i << " → " << j << "):\n";

                // outgoing arcs
                for (int id = 0; id<edges.size(); id++) {
                    if (edges[id].first == k) {
                        auto it = m_var_f_e_.find({id, {i, j}});
                        if (it != m_var_f_e_.end() && it->second != nullptr) {
                            constraint->SetCoefficient(it->second, 1);
                        } else {
                            std::cerr << "Warning: Variable for arc " << id << " and demand (" << i << ", " << j << ") not found.\n";
                        }

                        //std::cout << "  +1 * f_" << id << "_(" << i << "," << j << ")\n";
                    }
                }

                // incoming arcs
                for (int id = 0; id<edges.size(); id++) {
                    if (edges[id].second == k) {
                        auto it = m_var_f_e_.find({id, {i, j}});
                        if (it != m_var_f_e_.end() && it->second != nullptr) {
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

void LPSolver::SetObjective()
{
    // === Objective: maximize alpha ===
    solver->MutableObjective()->SetCoefficient(alpha, 1);
    solver->MutableObjective()->SetMinimization();
}


double LPSolver::getMaximumCongestion(const IGraph& graph) const {
    double max = 0;
    // === Print the solution ===
    std::unordered_map<int, double> total_flow_per_arc;

    for (const auto& [key, var] : m_var_f_e_) {
        if (!var) continue;
        double val = var->solution_value();
        int arc_id = std::get<0>(key);
        total_flow_per_arc[arc_id] += val;
    }

    for (int id = 0; id<edges.size(); id++) {
        int u = edges[id].first, v = edges[id].second;
        int rev_edge = -1;
        double total_flow = total_flow_per_arc[id];

        // to ensure that flow along anti-parallel arcs
        // will be added as absolute flow to its corresponding undirected link
        if (u > v) {
            auto rev_it = std::find(edges.begin(), edges.end(), std::make_pair(v, u));
            if(rev_it != edges.end()) {
                rev_edge = std::distance(edges.begin(), rev_it);
                total_flow = total_flow_per_arc[rev_edge];
                v = edges[id].first;
                u = edges[id].second;
            }
        }

        double capacity = graph.getEdgeCapacity(u, v);

        if(max < total_flow/capacity) {
            max = total_flow/capacity;
        }

    }


    return max;
}


void LPSolver::GetRoutingTable(const IGraph& graph) {
    std::cout << "\n=== Oblivious Routing Table ===\n";

    for (const auto& [key, var] : m_var_f_e_) {

        if (!var) {
            std::cerr << "Warning: variable pointer is null for key.\n";
            continue;
        }
        double value = var->solution_value();
        const auto& arc = std::get<0>(key);
        const std::pair<int, int>& demand = std::get<1>(key);




    }
}

void LPSolver::PrintSolution(const IGraph &graph) {

    double max_cong(0);
    // === Print the solution ===
    std::unordered_map<int, double> total_flow_per_arc;

    for (const auto& [key, var] : m_var_f_e_) {
        if (!var) continue;
        double val = var->solution_value();
        int arc_id = std::get<0>(key);
        total_flow_per_arc[arc_id] += val;
    }

    for (int id = 0; id<edges.size(); id++) {
        int u = edges[id].first, v = edges[id].second;

        if (u > v) continue;

        int rev_edge = -1;
        double total_flow = total_flow_per_arc[id];

        // to ensure that flow along anti-parallel arcs
        // will be added as absolute flow to its corresponding undirected link
        auto rev_it = std::find(edges.begin(), edges.end(), std::make_pair(v, u));
        if(rev_it != edges.end()) {
            rev_edge = std::distance(edges.begin(), rev_it);
            total_flow += total_flow_per_arc[rev_edge];
        }

        double capacity = graph.getEdgeCapacity(u, v);

        if(max_cong < total_flow/capacity) {
            max_cong = total_flow/capacity;
        }

    }
    std::cout << "Max congestion across all undirected arcs: " << max_cong << std::endl;

    PrintCommoditiesPerEdge(graph);
}

void LPSolver::PrintCommoditiesPerEdge(const IGraph& graph) {
    std::cout << "\n=== Commodities per Edge ===\n";

    // Iterate over all edges
    for (int edge_id = 0; edge_id < edges.size(); edge_id++) {
        const auto& edge = edges[edge_id];
        std::cout << "Edge " << edge_id << " (" << edge.first << " -> " << edge.second << "):\n";

        // Iterate over all demands for this edge
        for (const auto& d : demands) {
            auto it = m_var_f_e_.find({edge_id, {d.first, d.second}});
            if (it != m_var_f_e_.end() && it->second != nullptr) {
                double flow_value = it->second->solution_value();
                if (flow_value > 1e-9) { // print only non-zero flows
                    std::cout << "  Commodity (" << d.first << " -> " << d.second
                              << "): " << flow_value << "\n";
                }
            }
        }
    }
}

void LPSolver::storeFlow(EfficientRoutingTable& table) {
    // store the flow for each commodity in f_st_e
    // f_st_e.resize(edges.size()/2); // note that in vector edges, edges are stored two directed arcs, thus for outputting we only need one direction

    for (int edge_id = 0; edge_id < edges.size(); edge_id++) {
        const auto& edge = edges[edge_id];

        // Iterate over all demands for this edge
        for (const auto& d : demands) {
            auto it = m_var_f_e_.find({edge_id, {d.first, d.second}});


            if (it != m_var_f_e_.end() && it->second != nullptr) {
                double flow_value = it->second->solution_value();
                if (flow_value > 1e-9) { // print only non-zero flows
                    int e = g->getEdgeId(edges[edge_id].first, edges[edge_id].second);
                    if(edge.first < edge.second) {

                        table.addFlow(e, d.first, d.second, flow_value);
                        //f_e_st[edge][{d.first, d.second}] =  flow_value;
                    } else {
                        //f_e_st[{edge.second, edge.first}][{d.first, d.second}] =  -flow_value;

                        int anti_edge = g->getAntiEdge(e);;
                        table.addFlow(anti_edge, d.first, d.second, -flow_value);
                    }
                }
            }
        }
    }
    if(debug) {
        // Print the stored flow for each edge
        for (int edge_id = 0; edge_id < edges.size(); edge_id++) {
            const auto& edge = edges[edge_id];
            std::cout << "Edge (" << edge.first << " -> " << edge.second << "):\n";

            for (const auto& d : demands) {
                double flow_value = table.getFlow(edge_id, d.first, d.second);
                if (std::abs(flow_value) > 1e-9) {
                    std::cout << "  Commodity (" << d.first << " -> " << d.second
                              << "): " << flow_value << "\n";
                }
            }
        }
    }
}

double LPSolver::getCongestion(DemandMap& _demands, IGraph& g) const{
    // compute for the given demands and the stored flow values in f_st_e the generated congestion
    double max_cong = 0.0;

    std::unordered_map<std::pair<int, int>, double, PairHash> total_edge_congestion;

    for(int id = 0; id<_demands.size(); id++) {
        auto [source, target] = _demands.getDemandPair(id);
        double value = _demands.getDemandValue(id);

        for (int id = 0; id < edges.size(); id++) {
            if (edges[id].first > edges[id].second) continue; // Skip anti-parallel arcs

            auto it = m_var_f_e_.find({id, {source, target}});
            if (it != m_var_f_e_.end() && it->second != nullptr) {
                double flow = it->second->solution_value();
                total_edge_congestion[{edges[id].first, edges[id].second}] += std::abs(flow)*value; // Sum the absolute flow values for each edge
            }
        }
    }

    for(auto& [edge, flow] : total_edge_congestion) {
        max_cong = std::max(max_cong, flow); // Update the maximum congestion
    }
    return max_cong;
    /*
    std::vector<double> total_edge_congestion(edges.size(), 0.0);
    double max_cong = 0.0;

    for(const auto [d, value] : _demands) {
        int source = d.first, target = d.second;

        for (int id = 0; id < edges.size(); id++) {
            if (edges[id].first > edges[id].second) continue; // Skip anti-parallel arcs

            auto it = f_e_st.find({id, {source, target}});
            if (it != f_e_st.end() && it->second != nullptr) {
                double flow = it->second->solution_value();
                total_edge_congestion[id] += std::abs(flow)*value; // Sum the absolute flow values for each edge
            }
        }
    }

    for(const auto& congestion : total_edge_congestion) {
        if (congestion > max_cong) {
            max_cong = congestion; // Update the maximum congestion
        }
    }
    return max_cong;*/
}

