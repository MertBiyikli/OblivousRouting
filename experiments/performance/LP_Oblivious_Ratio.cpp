//
// Created by Mert Biyikli on 15.08.25.
//

#include "LP_Oblivious_Ratio.h"
#include "../../src/utils/hash.h"
#include <fstream>

void ObliviousRatio::init(Graph &g, const std::unordered_map<std::pair<int,int>, std::unordered_map<std::pair<int,int>, double, PairHash>, PairHash>& routing) {

    this->g = g;
    this->routing = routing;
    // Collect unique commodity pairs
    for (const auto& [edge, flow_map] : routing) {
        for (const auto& [st, _] : flow_map) {
            demands.insert(st);
        }
    }
}

double ObliviousRatio::solve() {
    double max_cong = -std::numeric_limits<double>::infinity();

    for (const auto& [edge, _] : routing) {
        auto result = solveWorstCaseDemandLPPerEdge(edge);
        max_cong = std::max(result, max_cong);
    }
    return max_cong;
}

double ObliviousRatio::solveWorstCaseDemandLPPerEdge(const std::pair<int, int> &edge) {
    std::unique_ptr<MPSolver> solver(MPSolver::CreateSolver("GLOP"));
    assert(solver);

    // Create variables
    std::unordered_map<std::pair<int, int>, MPVariable*,  PairHash> d_st;
    for (const auto& st : demands) {
        d_st[st] = solver->MakeNumVar(0.0, solver->infinity(), "d_" + std::to_string(st.first) + "_" + std::to_string(st.second));
    }


    MPVariable* alpha = solver->MakeNumVar(0.0, solver->infinity(), "z");

    auto und_edge = (edge.first< edge.second? edge : std::make_pair(edge.second, edge.first));
    double capacity = g.getEdgeCapacity(und_edge.first, und_edge.second);  // ensure you use the undirected form

    MPConstraint* cons = solver->MakeRowConstraint( 0,0.0,
        "cong_" + std::to_string(edge.first) + "_" + std::to_string(edge.second));

    for (const auto& [st, coeff] : routing.at(edge)) {
        cons->SetCoefficient(d_st[st], std::abs(coeff) / capacity);
    }

    cons->SetCoefficient(alpha, -1.0);

    for (const auto& [st, coeff] : routing.at(edge)) {
        MPConstraint* cons = solver->MakeRowConstraint( 0,solver->infinity(),
        "d_" + std::to_string(st.first) + "_" + std::to_string(st.second));

        cons->SetCoefficient(d_st[st], std::abs(coeff));
        cons->SetUB(capacity);
    }



    // Maximize z
    MPObjective* obj = solver->MutableObjective();
    obj->SetCoefficient(alpha, 1.0);
    obj->SetMaximization();

    const MPSolver::ResultStatus status = solver->Solve();

    if (status != MPSolver::OPTIMAL) {
        std::cerr << "Failed to solve LP\n";
        // print error status
        std::cerr << "Status: " << status << "\n";

        return -1.0;
    }else{
        return alpha->solution_value();
    }
}

