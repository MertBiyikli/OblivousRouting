//
// Created by Mert Biyikli on 30.09.25.
//

#include "electrical.h"
#include "../utils/my_math.h"
#include <random>
/*
    Main algorithm implementation
 */

void ElectricalMWU::init( bool debug,  boost::property_tree::ptree _params)
{

    n = graph.getNumNodes();
    m = graph.getNumEdges()/2;

    // set algorithm parameters
    roh = std::sqrt(2.0*static_cast<double>(m)); // Initialize roh based on the number of nodes
    alpha_local = std::log2(n)*std::log2(n); // Initialize alpha_local based on the number of nodes
    this->cap_X = m;
    this->iteration_count = 8.0*roh*std::log(m)/alpha_local;
    this->inv_m = 1.0 / static_cast<double>(m);
    this->x_fixed = 0;


    initEdgeDistances();
    buildIncidence();

    // init AMG
    amg = std::make_unique<AMGSolver>();
    if (amg == nullptr) {
        std::cerr << "Failed to create AMG solver instance.\n";
        return;
    }

    if (_params.empty()) {
        // set default parameters if not provided
        boost::property_tree::ptree config;
        boost::property_tree::read_json("../../configs/amg_configs.json", config);
        _params = config;
    }
    // tor parsing the configuration file for the AMG solver, e.g. coarsening and relaxation types
    amg->setSolverParams(_params);
    amg->init(graph, edge_weights, n, edges, debug);


    // compute Sketch matrix
    SketchMatrix = getSketchMatrix(m, n, 0.5);
    SketchMatrix_t = SketchMatrix.transpose();

    Eigen::MatrixXd UCt = SketchMatrix_t; // m × ℓ
    for (int e = 0; e < m; ++e)
        UCt.row(e) *= edge_capacities[e];

    auto X_dense = (B.transpose() * UCt); // n × ℓ
    X = X_dense.sparseView();

}

void ElectricalMWU::initEdgeDistances() {
    extractEdges();

    // note that the edges are stored undirected
    // --- Per-edge init ---
    edge_distances.assign(m, 1.0);
    edge_capacities.resize(m);
    edge_probabilities.resize(m);
    edge_weights.resize(m);

    for (int e = 0; e < m; ++e) {
        double cap = graph.getEdgeCapacity(e);      // undirected capacity accessor
        edge_capacities[e]  = cap;
        edge_probabilities[e] = edge_distances[e] / cap_X;
        edge_weights[e] = std::pow(cap, 2) / (edge_probabilities[e] + inv_m);
    }


}


void ElectricalMWU::extractEdges() {
    // in edges we only store the undirected edges
    edges.reserve(m);
    for (int e = 0; e < graph.getNumEdges(); e++) {
        auto [u, v] = graph.edgeEndpoints(e);
        if (u < v) {
            edges.emplace_back(u,v);
        }
    }

    // sort the edges based on the first node, then second node
    std::sort(edges.begin(), edges.end(),
              [](const std::pair<int,int> &a, const std::pair<int,int> &b) {
                  if (a.first != b.first) return a.first < b.first;
                  return a.second < b.second;
              });
}



void ElectricalMWU::buildIncidence()
{
    if (B.nonZeros() == 0) {
        B = Eigen::SparseMatrix<double>(m, n);
        std::vector<Eigen::Triplet<double>> T; T.reserve(2*m);
        for (int e = 0; e < m; ++e) {
            auto [u,v] = edges[e]; // u < v
            T.emplace_back(e, u, -1.0);
            T.emplace_back(e, v, +1.0);
        }
        B.setFromTriplets(T.begin(), T.end());
    }
}


void ElectricalMWU::updateEdgeDistances(const std::vector<double>& load) {
    for (int e = 0; e < m; ++e) {
        if (load[e] <= 0.0 || std::isnan(load[e])) continue;
        edge_distances[e] *= (1 + (1/(2*roh)) * load[e]);
        cap_X += edge_distances[e];
    }

    for (int e = 0; e < m; ++e) {
        double x = edge_distances[e];
        if (x <= 0.0 || std::isnan(x)) continue;

        double p = x / cap_X;
        edge_probabilities[e] = p;

        double cap = edge_capacities[e];
        double w   = std::pow(cap, 2) / (p + inv_m);
        edge_weights[e] = w;
    }


    amg->updateAllEdges(edge_weights, edges);
    amg->updateSolver();

}




void ElectricalMWU::run(LinearRoutingTable &table) {

    Eigen::VectorXd rhs = Eigen::VectorXd::Zero(n);
    Eigen::VectorXd potentials(n);
    std::vector<double> load(m, 0.0);

    auto start_run = timeNow();

    for (int t = 0; t < this->iteration_count; ++t) {
        auto start_oracle = timeNow();
        // --- main loop over sources (u -> x_fixed) ---
        for (int u = 0; u < n; ++u) {
            if (u == x_fixed) continue;

            // build RHS (sum is zero automatically)
            rhs.setZero();
            rhs[u]       =  1.0;
            rhs[x_fixed] = -1.0;

            auto start_oracle_pure = timeNow();
            potentials = amg->solve(rhs, epsilon_L);
            pure_oracles_running_times.push_back(duration(timeNow()-start_oracle_pure));
            addFlowToTable(u, potentials, table);
        }

        getApproxLoad(load);
        updateEdgeDistances(load); // keeps 'weights' vector in sync

        oracle_running_times.push_back(duration(timeNow() - start_oracle));
    }

    this->solve_time = duration(timeNow() - start_run);
}


void ElectricalMWU::addFlowToTable(const int& u, Eigen::VectorXd& potential, LinearRoutingTable &table) {
    // f_e(u) = w_e * (x[u]-x[v])
    auto t0 = timeNow();
    for (int e = 0; e < m; ++e) {
        double fval = edge_weights[e] * (potential[edges[e].first] - potential[edges[e].second]);

        if (std::abs(fval) > EPS) {
            // if value is negative, we push along the anti-edge
            const auto& [head, tail] = edges[e];

            // note that the edges in this loop iterates over all undirected edges, whereas we store the negative
            // as the anti-edge in the graph, so we need to get the original edge id and then get the anti-edge if needed
            const int& original_edge_id = graph.getEdgeId(head, tail);

            int e = (fval < 0 ? graph.getAntiEdge(original_edge_id) : original_edge_id);
            t0 = timeNow();
            table.addFlow(e, u, std::abs(fval));
            this->transformation_time += duration(timeNow() - t0);
        }
    }
}


void ElectricalMWU::scaleFlowDown(LinearRoutingTable& table) {
    // scale the flow from the adjacency list flow
    if (iteration_count > 0) {
        const double inv_iters = 1.0 / static_cast<double>(iteration_count);
        for (int e = 0; e < graph.getNumEdges(); ++e) // dont use m here. m is undirected edges only
            for (double &val : table.src_flows[e]) val *= (inv_iters);
    }
}

void ElectricalMWU::getApproxLoad(std::vector<double>& load) {
    const int ell = X.cols();
    Eigen::VectorXd rhs(n), sol(n);
    Eigen::VectorXd d(m);

    edge_diffs.resize(static_cast<size_t>(m) * ell);
    double K_obs = 0.0;

    for (int i = 0; i < ell; ++i) {
        rhs = X.col(i);
        sol = amg->solve(rhs, epsilon_L);

        for (int e = 0; e < m; ++e) {
            const auto [u,v] = edges[e];        // B has -1 at u, +1 at v
            d[e] = sol[v] - sol[u];             // signed diff consistent with B
            edge_diffs[size_t(e) * ell + i] = std::abs(d[e]); // keep abs for median
        }
    }

    // update the Laplacian error bound K based on the observed max diff in the sketch space
    if (std::abs(K-1.0) < EPS) {
            Eigen::VectorXd y =SketchMatrix_t * d;
            K_obs = std::max(K_obs, y.cwiseAbs().maxCoeff());
    }
    K = std::max(K, 1.5 * K_obs); // safety factor, monotone

    if (K > 0.0) {
        epsilon_L = epsilon / (8.0 * m * std::pow(n, 4) * K);
        epsilon_L = std::max(epsilon_L, 1e-12);
    }

    // recover norm
    for (int e = 0; e < m; ++e) {
        double* __restrict arr = &edge_diffs[static_cast<size_t>(e) * ell];
        std::nth_element(arr, arr + (ell >> 1), arr + ell);
        const double med = arr[ell >> 1];
        load[e] = edge_weights[e] * med;
    }
}



Eigen::MatrixXd ElectricalMWU::getSketchMatrix(int _m, int _n, double eps) {
    double c = 1.1;
    double delta = NegativeExponent(_n, 10); // Set delta to a small value or a default value
    double epsilon = eps; // Set epsilon to the provided value or a default value

    int l = (c/(epsilon*epsilon) * std::log(1.0/delta) );

    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::cauchy_distribution<double> dist(0.0, 1.0);

    Eigen::MatrixXd C(l, _m);
    double val = 0;
    for (int i = 0; i < l; ++i) {
        for (int j = 0; j < _m; ++j) {
            C.coeffRef(i, j) = val;
        }
    }

    return C;
}