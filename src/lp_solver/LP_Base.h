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
#include "../datastructures/GraphCSR.h"


using namespace operations_research;

class LP : public ObliviousRoutingSolver{
public:
    /*
    const GraphCSR *g = nullptr;
    bool debug = false;
    int n, m;
    std::unique_ptr<MPSolver> solver;
    MPVariable* alpha;
    std::vector<Demand> demands;
    std::vector<std::pair<int, int> > edges;

    LP() : solver(nullptr), alpha(nullptr) {}
    virtual ~LP() = default;

    virtual void CreateVariables(const GraphCSR& graph) = 0;
    virtual void CreateConstraints(const GraphCSR& graph) = 0;
    virtual void SetObjective() = 0;
    virtual void PrintSolution(const GraphCSR& graph) = 0;
    void init(const GraphADJ& graph) {
        edges.clear();
        n = graph.getNumNodes();
        m = graph.getNumEdges()*2; // since we store each undirected edge as two anti parallel directed edges

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

    void init(const GraphCSR& graph) {
        edges.clear();
        n = graph.getNumNodes();
        m = graph.getNumEdges(); // *2 , Not here

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

    void runSolve(const IGraph &graph_) override {
        // Safe downcast: IGraph -> Graph_csr
        auto* csr = dynamic_cast<const GraphCSR*>(&graph_);
        if (!csr) {
            throw std::runtime_error("LP::runSolve: expected Graph_csr as input graph.");
        }

        g = csr;  // non-owning pointer, just an alias

        if (Run(*csr)) {
            if (debug) {
                this->PrintSolution(*csr);
            }
        } else {
            std::cout << "ERROR: Failed to complete LP." << std::endl;
        }
    }


    bool Run(const GraphCSR& graph) {

        init(graph);
        CreateVariables(graph);
        CreateConstraints(graph);
        SetObjective();
        // === Solve the LP ===
        if (solver->Solve() != MPSolver::OPTIMAL) {
            std::cerr << "WARNING: LP is infeasible or unbounded.\n";
            return false;
        } else {
            std::cout << "Optimal offline congestion: " << alpha->solution_value() << std::endl;

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
    */
};

#endif //OBLIVIOUSROUTING_LP_BASE_H
