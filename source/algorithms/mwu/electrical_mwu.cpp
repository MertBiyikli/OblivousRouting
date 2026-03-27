//
// Created by Mert Biyikli on 20.03.26.
//

#include "../../../include/algorithms/mwu/electrical_mwu.h"
#include "../../../include/utils/my_math.h"
#include <random>

boost::property_tree::ptree make_amg_params() {
    boost::property_tree::ptree pt;

    // ---- solver ----
    pt.put("solver.type", "cg");
    pt.put("solver.tol", 1e-8);
    pt.put("solver.maxiter", 100);

    // AMG details (typical)
    pt.put("precond.coarsening.type", "smoothed_aggr_emin");
    pt.put("precond.relax.type", "spai1");

    return pt;
}

void ElectricalMWU::init( bool debug,  boost::property_tree::ptree _params)
{
    auto t0 = timeNow();
    n = graph.getNumNodes();
    m = graph.getNumUndirectedEdges();

    // set algorithm parameters
    roh = std::sqrt(2.0*static_cast<double>(m));
    alpha_local = std::log2(n)*std::log2(n);
    this->cap_X = m;
    this->iteration_count = std::max(1, (int)std::ceil(8.0 * roh * std::log((double)m) / alpha_local));
    this->inv_m = 1.0 / static_cast<double>(m);
    this->x_fixed = 0;

    initEdgeDistances();

    boost::property_tree::ptree params = make_amg_params();
    initAMGSolver(params);

    // compute Sketch matrix
    SketchMatrix_t = getSketchMatrix(m, n, 0.5).transpose();

    Eigen::MatrixXd UCt = SketchMatrix_t; // m × ℓ
    for (int e = 0; e < m; ++e)
        UCt.row(e) *= edge_capacities[e];

    auto B = buildIncidence();
    X = (B.transpose() * UCt).sparseView(); // n × ℓ

    solve_time += duration(timeNow()-t0);
}


void ElectricalMWU::initAMGSolver(boost::property_tree::ptree _params) {
    // init AMG
    amg = std::make_unique<LaplacianSolver>();
    if (amg == nullptr) {
        std::cerr << "Failed to create AMG solver instance.\n";
        return;
    }
    // tor parsing the configuration file for the AMG solver, e.g. coarsening and relaxation types
    amg->setSolverParams(_params);
    amg->init(graph, edge_weights, n, edges, debug);
}

/*
 * The main loop of the MWU algorithm.
 * For each iteration, we loop over all sources (except the fixed one), build the RHS for the Laplacian system,
 * solve for potentials, and then add the computed flow to the routing table.
 * After processing all sources, we compute the approximate load and update edge distances accordingly.
 */
void ElectricalMWU::run(LinearRoutingTable &table) {

    auto t0 = timeNow();
    Eigen::VectorXd rhs = Eigen::VectorXd::Zero(n);
    Eigen::VectorXd potentials(n);
    std::vector<double> load(m, 0.0);
    solve_time += duration(timeNow() - t0);


    for (int t = 0; t < this->iteration_count; ++t) {

        auto oracle_iteration = 0.0;
        // --- main loop over sources (u -> x_fixed) ---
        for (int u = 0; u < n; ++u) {
            if (u == x_fixed) continue;

            t0 = timeNow();
            rhs.setZero();
            rhs[u]       =  1.0;
            rhs[x_fixed] = -1.0;
            double setup_time = duration(timeNow() - t0);

            auto start_oracle_pure = timeNow();
            potentials = amg->solve(rhs, epsilon_L);
            double oracle_time_iter = duration(timeNow() - start_oracle_pure);
            oracle_iteration += oracle_time_iter;

            // addFlowToTable measures its own time and adds to transformation_time
            addFlowToTable(u, potentials, table);

            // solve_time includes setup_time + oracle_time (but not transformation_time)
            solve_time += setup_time + oracle_time_iter;
        }

        t0 = timeNow();
        getApproxLoad(load);
        updateDistances(load);
        double weight_update_time = duration(timeNow() - t0);
        mwu_weight_update_time += weight_update_time;
        solve_time += weight_update_time;

        oracle_running_times.push_back(oracle_iteration);
    }


}


/*
 * This function takes the potentials obtained from solving the Laplacian system
 * and computes the flow on each edge based on the potential difference and edge resistances.
 */
void ElectricalMWU::addFlowToTable(const int& u, Eigen::VectorXd& potential, LinearRoutingTable &table) {
    auto t0 = timeNow();
    for (int e = 0; e < m; ++e) {
        double fval = edge_weights[e] * (potential[edges[e].first] - potential[edges[e].second]);

        if (std::abs(fval) > SOFT_EPS) {
            // if value is negative, we push along the anti-edge
            const auto& [head, tail] = edges[e];

            // note that the edges in this loop iterates over all undirected edges, whereas we store the negative
            // as the anti-edge in the graph, so we need to get the original edge id and then get the anti-edge if needed
            const int& original_edge_id = graph.getEdgeId(head, tail);

            int e_orig = (fval < 0 ? graph.getAntiEdge(original_edge_id) : original_edge_id);

            table.addFlow(e_orig, u, std::abs(fval));
        }
    }
    this->transformation_time += duration(timeNow() - t0);
}


/*
 * This function computes the approximate load on each edge based on the current potentials obtained from solving the Laplacian system.
 * It uses the sketching matrix to project the flow differences into a lower-dimensional space and updates the approximate load estimates accordingly.
 */
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
            const auto& [u,v] = edges[e];
            d[e] = sol[v] - sol[u];             // signed diff consistent with B
            edge_diffs[size_t(e) * ell + i] = std::abs(d[e]); // keep abs for median
        }

        if (!K_initialized) {
            Eigen::VectorXd y = SketchMatrix_t.transpose() * d; // (ℓ×m)*(m) = ℓ   (clearer than using transpose)
            K_obs = std::max(K_obs, y.cwiseAbs().maxCoeff());
        }
    }


    if (!K_initialized) {
        K = std::max(K, 1.5 * K_obs); // safety factor, monotone
        K_initialized = true;
    }

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


/*
 * This function updates the edge distances based on the computed load. It iterates over all edges,
 * and for those with positive load, it updates the edge distance using the formula: x_e *= (1 + (1/(2*roh)) * load[e]).
 */
void ElectricalMWU::updateDistances(const std::vector<double>& load) {
    cap_X = 0.0;
    for (int e = 0; e < m; ++e) {
        const double le = load[e];
        if (le > 0.0 && std::isfinite(le)) {
            edge_distances[e] *= (1.0 + (1.0/(2.0*roh)) * le);
        }
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

/*
 * After computing the flows in the main loop, we need to scale them down by the number of iterations to get the average flow.
 */
void ElectricalMWU::scaleFlowDown(LinearRoutingTable& table) {
    // scale the flow from the adjacency list flow
    auto start_transfo = timeNow();
    if (iteration_count > 0) {
        const double inv_iters = 1.0 / static_cast<double>(iteration_count);
        for (int e = 0; e < graph.getNumDirectedEdges(); ++e) // dont use m here. m is undirected edges only
            for (double &val : table.src_flows[e]) val *= (inv_iters);
    }
    this->transformation_time += duration(timeNow() - start_transfo);
}

void ElectricalMWU::initEdgeDistances() {
    extractEdges();

    // note that the edges are stored undirected-
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
    for (int e = 0; e < graph.getNumDirectedEdges(); e++) {
        auto [u, v] = graph.getEdgeEndpoints(e);
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

Eigen::SparseMatrix<double> ElectricalMWU::buildIncidence()
{
    Eigen::SparseMatrix<double> B(m, n);
    std::vector<Eigen::Triplet<double>> T; T.reserve(2*m);
    for (int e = 0; e < m; ++e) {
        auto [u,v] = edges[e]; // u < v
        T.emplace_back(e, u, -1.0);
        T.emplace_back(e, v, +1.0);
    }
    B.setFromTriplets(T.begin(), T.end());
    return B;
}


/*
 * Sketching matrix generation using Cauchy distribution for L1 norm approximation.
 * The goal is to create a matrix that can be used to project the flow differences into a lower-dimensional space while preserving the L1 norm properties.
 */
Eigen::MatrixXd ElectricalMWU::getSketchMatrix(int _m, int _n, double eps) {
    double c = 1.1;
    double delta = NegativeExponent(_n, 10);
    double epsilon = eps;

    int l = (c/(epsilon*epsilon) * std::log(1.0/delta) );

    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::cauchy_distribution<double> dist(0.0, 1.0);

    Eigen::MatrixXd C(l, _m);
    double val = 0;
    for (int i = 0; i < l; ++i) {
        for (int j = 0; j < _m; ++j) {
            val = dist(gen);
            C.coeffRef(i, j) = val;
        }
    }
    return C;
}