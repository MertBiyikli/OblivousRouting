//
// Created by Mert Biyikli on 20.03.26.
//

#include "algorithms/oblivious/mwu/electrical_mwu.h"
#include "utils/my_math.h"
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

Result<void> ElectricalMWU::init( bool debug,  boost::property_tree::ptree _params)
{
    auto t0 = timeNow();
    n = graph.getNumNodes();
    m = graph.getNumUndirectedEdges();

    // set algorithm parameters
    roh = std::sqrt(2.0*static_cast<double>(m));
    alpha_local = std::log2(n)*std::log2(n);
    this->cap_X = m;
    metrics.iteration_count = std::max(1, (int)std::ceil(8.0 * roh * std::log((double)m) / alpha_local));
    this->inv_m = 1.0 / static_cast<double>(m);
    this->x_fixed = 0;

    initEdgeDistances();

    boost::property_tree::ptree params = make_amg_params();
    auto amg_init = initAMGSolver(params);
    if (!amg_init) {
        return getError(amg_init);
    }

    if ( use_sketching ) {
        // compute Sketch matrix
        SketchMatrix_t = getSketchMatrix(0.5).transpose();

        Eigen::MatrixXd UCt = SketchMatrix_t; // m × ℓ
        for (int e = 0; e < m; ++e)
            UCt.row(e) *= edge_capacities[e];

        auto B = buildIncidence();
        X = (B.transpose() * UCt).sparseView(); // n × ℓ
    }
    metrics.solve_time += duration(timeNow()-t0);
    return {};
}


Result<void> ElectricalMWU::initAMGSolver(boost::property_tree::ptree _params) {
    // init AMG
    amg = std::make_unique<LaplacianSolver>();
    if (!amg) {
       return makeErrorMessage(ErrorCode::InvalidSolver, "Failed to create AMG solver instance.");
    }
    // tor parsing the configuration file for the AMG solver, e.g. coarsening and relaxation types
    amg->setSolverParams(_params);
    amg->init(graph, edge_weights, n, edges, debug);
    return {};
}

/*
 * The main loop of the MWU algorithm.
 * For each iteration, we loop over all sources (except the fixed one), build the RHS for the Laplacian system,
 * solve for potentials, and then add the computed flow to the routing table.
 * After processing all sources, we compute the approximate load and update edge distances accordingly.
 */
Result<void> ElectricalMWU::run(LinearRoutingTable &table) {

    Eigen::VectorXd rhs = Eigen::VectorXd::Zero(n);
    std::vector<double> load(m, 0.0);

    Result<Eigen::VectorXd> potentials(n);


    for (int t = 0; t < metrics.iteration_count; ++t) {

        auto t0 = timeNow();
        auto oracle_iteration = 0.0;
        // --- main loop over sources (u -> x_fixed) ---
        for (int u = 0; u < n; ++u) {
            if (u == x_fixed) continue;

            rhs.setZero();
            rhs[u]       =  1.0;
            rhs[x_fixed] = -1.0;

            potentials = amg->solve(rhs, epsilon_L);

            if (!potentials) {
                return getError(potentials);
            }
            double oracle_time_iter = duration(timeNow() - t0);
            oracle_iteration += oracle_time_iter;

            // addFlowToTable measures its own time and adds to transformation_time
            addFlowToTable(u, potentials.value(), table);

            // solve_time includes setup_time + oracle_time (but not transformation_time)
            metrics.solve_time += oracle_time_iter;
        }

        t0 = timeNow();

        // Compute the load on the edges
        Result<void> load_comp;
        if (use_sketching) {
            load_comp = getApproxLoad(load);
        }else {
            load_comp = getExactLoad(load);
        }

        if (!load_comp) {
            return getError(load_comp);
        }
        metrics.load_computation_time += duration(timeNow()-t0);

        t0 = timeNow();
        auto update = updateDistances(load);
        if (!update) {
            return getError(update);
        }
        double weight_update_time = duration(timeNow() - t0);
        metrics.mwu_weight_update_time += weight_update_time;

        metrics.oracle_running_times.push_back(oracle_iteration);
    }
    return {};
}

/*
 * This function takes the potentials obtained from solving the Laplacian system
 * and computes the flow on each edge based on the potential difference and edge resistances.
 */
void ElectricalMWU::addFlowToTable(const int& source,
                                   const Eigen::VectorXd& potential,
                                   LinearRoutingTable& table) {
    auto t0 = timeNow();

    for (int e = 0; e < m; ++e) {
        const auto [a, b] = edges[e];

        /*
         * Important convention:
         *
         * Your incidence matrix uses B[e,a] = -1 and B[e,b] = +1.
         * Therefore B * phi on edge e is phi[b] - phi[a].
         *
         * If signed_flow > 0, the electrical current corresponds to
         * source-side flow from b to a in the routing table.
         *
         * If signed_flow < 0, it corresponds to flow from a to b.
         */
        const double signed_flow =
            edge_weights[e] * (potential[b] - potential[a]);

        if (!std::isfinite(signed_flow)) {
            continue;
        }

        if (std::abs(signed_flow) <= EPS) {
            continue;
        }

        int from;
        int to;
        const double amount = std::abs(signed_flow);

        if (signed_flow > 0.0) {
            // Positive B-flow: store b -> a.
            from = b;
            to   = a;
        } else {
            // Negative B-flow: store a -> b.
            from = a;
            to   = b;
        }

        const int directed_edge_id = graph.getEdgeId(from, to);


        table.addFlow(directed_edge_id, source, amount);

    }

    metrics.transformation_time += duration(timeNow() - t0);
}

void ElectricalMWU::setEpsilon(double eps) {
    this->epsilon_L = eps;
}


/*
 * This function computes the approximate load on each edge based on the current potentials obtained from solving the Laplacian system.
 * It uses the sketching matrix to project the flow differences into a lower-dimensional space and updates the approximate load estimates accordingly.
 */
Result<void> ElectricalMWU::getApproxLoad(std::vector<double>& load) {
    const int ell = X.cols();
    Eigen::VectorXd rhs(n);
    Eigen::VectorXd d(m);

    edge_diffs.resize(static_cast<size_t>(m) * ell);
    double K_obs = 0.0;

    for (int i = 0; i < ell; ++i) {
        rhs = X.col(i);
        auto sol = amg->solve(rhs);

        if (!sol) {
            return getError(sol);
        }

        for (int e = 0; e < m; ++e) {
            const auto& [u,v] = edges[e];
            d[e] = sol.value()[v] - sol.value()[u];             // signed diff consistent with B
            edge_diffs[size_t(e) * ell + i] = std::abs(d[e]); // keep abs for median
        }
    }

    if (!K_initialized) {
        Eigen::VectorXd y = SketchMatrix_t.transpose() * d; // (ℓ×m)*(m) = ℓ   (clearer than using transpose)
        K_obs = std::max(K_obs, y.cwiseAbs().maxCoeff());
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
    return {};
}

/**
 * Compute exact loads without sketching approximation.
 * For each demand edge f, we solve the Laplacian with unit current injection/extraction at the endpoints.
 * Then for each edge e, we accumulate the load contribution: load_w(e) += w_e * |b_e L† b_f^T|
 * where b_e L† b_f^T is the potential difference across edge e when unit current is injected on edge f.
 * Final formula: load_w(e) = w_e * Σ_f |b_e L† b_f^T|
 */
Result<void> ElectricalMWU::getExactLoad(std::vector<double>& load) {
    Eigen::VectorXd rhs(n);
    // Initialize load to zero
    std::fill(load.begin(), load.end(), 0.0);

    // For each demand edge f
    for (int f = 0; f < m; ++f) {
        // Create RHS for unit current demand on edge f
        // Unit current source at u_f and sink at v_f
        rhs.setZero();
        const auto& [u_f, v_f] = edges[f];
        rhs[u_f] = 1.0;
        rhs[v_f] = -1.0;

        // Solve Laplacian system: L * potential = rhs
        auto sol = amg->solve(rhs);
        if (!sol) {
            return getError(sol);
        }

        // For each edge e, compute load contribution from this demand
        // load contribution = w_e * |potential_diff_e|
        // where potential_diff_e = b_e^T * L† * b_f = sol[v_e] - sol[u_e]
        for (int e = 0; e < m; ++e) {
            const auto& [u_e, v_e] = edges[e];
            double pot_diff = sol.value()[v_e] - sol.value()[u_e];  // b_e^T L† b_f^T
            load[e] += edge_weights[e] * std::abs(pot_diff);
        }
    }
    return {};
}

/*
 * This function updates the edge distances based on the computed load. It iterates over all edges,
 * and for those with positive load, it updates the edge distance using the formula: x_e *= (1 + (1/(2*roh)) * load[e]).
 */
Result<void> ElectricalMWU::updateDistances(const std::vector<double>& load) {
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


    auto update = amg->updateAllEdges(edge_weights, edges);
    if (!update) {
        return getError(update);
    }
    amg->updateSolver();

    return {};
}

/*
 * After computing the flows in the main loop, we need to scale them down by the number of iterations to get the average flow.
 */
Result<void> ElectricalMWU::scaleFlowDown(LinearRoutingTable& table) {
    // scale the flow from the adjacency list flow
    auto time_transfo = timeNow();

    if (metrics.iteration_count > 0) {
        const double inv_iters = 1.0 / static_cast<double>(metrics.iteration_count);
        for (int e = 0; e < graph.getNumDirectedEdges(); ++e) // dont use m here. m is undirected edges only
            for (double &val : table.src_flows[e]) val *= (inv_iters);
    }else {
        return makeErrorMessage(ErrorCode::NumericalFailure, "Dividing by zero.");
    }
    metrics.transformation_time += duration(timeNow() - time_transfo);
    return {};
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
    std::vector<Eigen::Triplet<double>> T;
    T.reserve(2*m);

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
Eigen::MatrixXd ElectricalMWU::getSketchMatrix(double eps) {
    double c = 1.1;
    double delta = NegativeExponent(n, 2);
    double epsilon = eps;

    int l = (c/(epsilon*epsilon) * std::log(1.0/delta) );

    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::cauchy_distribution<double> dist(0.0, 1.0);

    Eigen::MatrixXd C(l, m);
    double val = 0;
    for (int i = 0; i < l; ++i) {
        for (int j = 0; j <m; ++j) {
            val = dist(gen);
            C.coeffRef(i, j) = val;
        }
    }
    return C;
}