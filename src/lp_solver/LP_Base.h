//
// Created by Mert Biyikli on 01.07.25.
//

#ifndef OBLIVIOUSROUTING_LP_BASE_H
#define OBLIVIOUSROUTING_LP_BASE_H

#include "../solver/solver.h"
#include "ortools/linear_solver/linear_solver.h"
#include <map>
#include <vector>
#include <tuple>
#include "../utils/hash.h"
#include "../solver/solver.h"
#include "../datastructures/IGraph.h"


using namespace operations_research;

class LP : public ObliviousRoutingSolver{
public:

    const IGraph* g = nullptr;
    bool debug = false;
    int n, m;
    std::unique_ptr<MPSolver> solver;
    MPVariable* alpha;
    std::vector<Demand> demands;
    std::vector<std::pair<int, int> > edges;

    LP() : solver(nullptr), alpha(nullptr) {}
    LP(IGraph& graph) {
        init(graph);
    }
    virtual ~LP() = default;

    virtual void CreateVariables(const IGraph& graph) = 0;
    virtual void CreateConstraints(const IGraph& graph) = 0;
    virtual void SetObjective() = 0;
    virtual void PrintSolution(const IGraph& graph) = 0;


    void init(const IGraph& graph) {
        edges.clear();
        n = graph.getNumNodes();
        m = graph.getNumDirectedEdges(); // *2 , Not here
        g = &graph;

        if (!solver) {
            solver.reset(MPSolver::CreateSolver("GLOP"));
        }
        if (!solver) {
            throw std::runtime_error("GLOP solver unavailable.");
        }

        for(int u = 0; u<graph.getNumNodes(); u++) {
            for(int v : graph.neighbors(u)) {
                edges.push_back({u, v});
            }
        }

        std::sort(edges.begin(), edges.end(), [](std::pair<int, int>& lhs, std::pair<int, int>& rhs) {
            if (lhs.first != rhs.first)
                return lhs.first < rhs.first;   // descending order for first
            return lhs.second < rhs.second;
        });
    }

    std::unique_ptr<RoutingScheme> solve() override {
        assert(g != nullptr);
        EfficientRoutingTable table;
        table.init(g->getNumDirectedEdges());

        // solver must fill table via computeBasisFlows
        if (Run(*g, table)) {
            if (debug) {
                this->PrintSolution(*g);
            }
        } else {
            std::cout << "ERROR: Failed to complete LP." << std::endl;
        }


        return std::make_unique<NonLinearRoutingScheme>(
      *g,
      std::move(table));
    }


    bool Run(const IGraph& graph, EfficientRoutingTable& table) {

        init(graph);
        CreateVariables(graph);
        CreateConstraints(graph);
        SetObjective();
        // === Solve the LP ===
        if (solver->Solve() != MPSolver::OPTIMAL) {
            std::cerr << "WARNING: LP is infeasible or unbounded.\n";
            return false;
        } else {
            // std::cout << "Optimal offline congestion: " << alpha->solution_value() << std::endl;
            storeFlow(table);
            return true;
        }
    }

    double getCongestion() {
        if (alpha) {
            return alpha->solution_value();
        } else {
            throw std::runtime_error("LP not solved yet, alpha variable is null.");
        }
    }

    void setDebug(bool debug) {this->debug = debug;}

    virtual void storeFlow(EfficientRoutingTable& table) = 0;

};

#endif //OBLIVIOUSROUTING_LP_BASE_H
