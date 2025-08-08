//
// Created by Mert Biyikli on 09.07.25.
//

#include "laplcian_solver.h"

LaplacianSolver::LaplacianSolver():
    colptr_julia(nullptr),
    rowval_julia(nullptr),
    nzval_julia(nullptr),
    bvec_julia(nullptr),
    arg1(nullptr),
    arg2(nullptr)
{
    initted = false;
}



void LaplacianSolver::init(const RaeckeGraph &_g, const std::unordered_map<std::pair<int, int>, double> &_w_edges2weights) {
    if(initted)
        return;

    m_graph = _g;
    m = _g.getNumEdges();
    n = _g.getNumNodes();

    for(const auto& [edge, weight] : _w_edges2weights) {
        w_edges2weights[edge] = weight;
    }

    if(debug) {
        if (jl_exception_occurred()) {
            jl_value_t *err = jl_exception_occurred();
            jl_function_t *showerror = jl_get_function(jl_base_module, "showerror");
            jl_value_t *msg = jl_call1(showerror, err);
            const char *errmsg = jl_string_ptr(msg);
            std::cerr << "Julia error: " << errmsg << std::endl;
            throw std::runtime_error("Julia load failure.");
        }
    }

    initted = true;

    // init bvec
    bvec.resize(n);
    bvec.setZero();

    if(debug) {
        std::cout << "Initializing LaplacianSolver with " << n << " nodes and " << m << " edges." << std::endl;
        std::cout << "Edge weights:\n";
        for (const auto &[edge, weight] : _w_edges2weights) {
            std::cout << "(" << edge.first << ", " << edge.second << ") : " << weight << std::endl;
        }
    }

    buildCSC();


    std::cout << "BEFORE memcpy: colptr.size() = " << colptr.size()
              << ", rowval.size() = " << rowval.size()
              << ", nzval.size() = " << nzval.size() << std::endl;


    // Prepare args
    jl_value_t *int_type = jl_apply_array_type((jl_value_t *) jl_int32_type, 1);
    jl_value_t *float_type = jl_apply_array_type((jl_value_t *) jl_float64_type, 1);

    // Allocate arrays in Julia and copy data
    colptr_julia = jl_alloc_array_1d(int_type, colptr.size());
    memcpy(jl_array_data(colptr_julia), colptr.data(), colptr.size() * sizeof(int32_t));

    rowval_julia = jl_alloc_array_1d(int_type, rowval.size());
    memcpy(jl_array_data(rowval_julia), rowval.data(), rowval.size() * sizeof(int32_t));

    nzval_julia = jl_alloc_array_1d(float_type, nzval.size());
    memcpy(jl_array_data(nzval_julia), nzval.data(), nzval.size() * sizeof(double));

    bvec_julia = jl_alloc_array_1d(float_type, bvec.size());
    memcpy(jl_array_data(bvec_julia), bvec.data(), bvec.size() * sizeof(double));


    std::cout << "n = " << n << std::endl;
    std::cout << "m = " << m << std::endl;
    std::cout << "colptr size = " << colptr.size() << std::endl;
    std::cout << "rowval size = " << rowval.size() << std::endl;
    std::cout << "nzval size = " << nzval.size() << std::endl;
    std::cout << "bvec size = " << bvec.size() << std::endl;

    // Prepare arguments for the Julia function
    arg1 = jl_box_int64(n);
    arg2 = jl_box_int64(m);

}



void LaplacianSolver::buildCSC() {
    // Prepare edge list
    std::vector<std::tuple<int, int, double>> edgeList;
    edgeList.reserve(w_edges2weights.size() * 2);

    for (const auto &[edge, weight] : w_edges2weights) {
        edgeList.emplace_back(edge.first, edge.second, weight);
    }

    // Sort edge list by column, then by row
    std::sort(edgeList.begin(), edgeList.end(),
              [](const auto &a, const auto &b) {
                  if (std::get<1>(a) != std::get<1>(b))
                      return std::get<1>(a) < std::get<1>(b);
                  return std::get<0>(a) < std::get<0>(b);
              });

    colptr.assign(n + 1, 0);
    rowval.clear();
    nzval.clear();

    int current_col = 0;
    int nnz_so_far = 0;

    for (const auto &[row, col, val] : edgeList) {
        while (col > current_col) {
            colptr[current_col + 1] = nnz_so_far;
            current_col++;
        }

        rowval.push_back(row);
        nzval.push_back(val);
        nnz_so_far++;
    }

    // Fill any trailing columns with current nnz
    while (current_col < n) {
        colptr[current_col + 1] = nnz_so_far;
        current_col++;
    }

    // Convert indices to 1-based for Julia
    for (auto &x : rowval) x += 1;
    for (auto &x : colptr) x += 1;

    if(debug) {
        std::cout << "LaplacianSolver: CSC representation built with " << colptr.size() << " columns, "
                  << rowval.size() << " rows, and " << nzval.size() << " non-zero values." << std::endl;

        // Debug prints
        std::cout << "colptr:";
        for (auto v: colptr) std::cout << " " << v;
        std::cout << std::endl;

        std::cout << "rowval:";
        for (auto v: rowval) std::cout << " " << v;
        std::cout << std::endl;

        std::cout << "nzval:";
        for (auto v: nzval) std::cout << " " << v;
        std::cout << std::endl;
    }
}



Eigen::VectorXd LaplacianSolver::solve(const Eigen::VectorXd &b) {

    // Build a dense matrix version for validation
    Eigen::MatrixXd Ldense = Eigen::MatrixXd::Zero(n, n);

    // Fill diagonals and off-diagonals correctly
    for (const auto& [edge, w] : w_edges2weights) {
        int u = edge.first;
        int v = edge.second;

        if (u == v) continue;  // ignore self-loops
        Ldense(u,u) += w;
        Ldense(v,v) += w;
        Ldense(u,v) -= w;
        Ldense(v,u) -= w;
    }

    // Check symmetry
    std::cout << "Symmetry check: " << ( (Ldense - Ldense.transpose()).norm() ) << std::endl;



    if(!initted) {
        throw std::runtime_error("LaplacianSolver not initialized. Call init() first.");
    }
    if(b.size() != n) {
        throw std::runtime_error("Input vector size does not match the number of nodes in the graph.");
    }

    // init the bvec if the dimensions are matching
    // TODO: try to identify this differently
    if(b.size() == m_graph.getNumNodes()) {
        bvec = b;
        // Copy bvec to Julia array
        memcpy(jl_array_data(bvec_julia), bvec.data(), bvec.size() * sizeof(double));
    } else {
        throw std::runtime_error("Input vector size does not match the number of nodes in the graph.");
    }
    jl_function_t* func = jl_get_function(jl_main_module, "laplacian_solve");
    if (func == nullptr) {
        throw std::runtime_error("Function not found.");
    }


    if (colptr.size() != static_cast<typeof(colptr.size())>(n + 1))
        throw std::runtime_error("colptr size must be n+1");

    if (rowval.size() != nzval.size())
        throw std::runtime_error("rowval and nzval sizes must match");

    if (!std::is_sorted(colptr.begin(), colptr.end()))
        throw std::runtime_error("colptr must be non-decreasing");


    // pass the debug flag as a Julia boolean
    jl_value_t* debug_julia = jl_box_bool(debug);

    jl_value_t* args[7] = {
            arg1,
            arg2,
            (jl_value_t *) colptr_julia,
            (jl_value_t *) rowval_julia,
            (jl_value_t *) nzval_julia,
            (jl_value_t *) bvec_julia,
            debug_julia
    };

    jl_array_t* result = (jl_array_t *)jl_call(func, args, 6);

    if(debug) {
        if (jl_exception_occurred()) {
            jl_value_t *e = jl_exception_occurred();
            if (e) {
                jl_function_t *sprint = jl_get_function(jl_base_module, "sprint");
                jl_function_t *showerror = jl_get_function(jl_base_module, "showerror");
                jl_value_t *str = jl_call2(sprint, showerror, e);
                std::cerr << "Julia error: " << jl_string_ptr(str) << std::endl;
            }
            throw std::runtime_error("Julia call failed.");
        }
    }
    // Convert the result to Eigen::VectorXd
    if(result) {
        double* data_ptr = (double*)jl_array_data(result);
        size_t len = jl_array_len(result);

        Eigen::VectorXd x(len);
        for (size_t i = 0; i < len; ++i) {
            x[i] = data_ptr[i];
        }

        return x;
    }else{
        if(debug) {
            jl_value_t *e = jl_exception_occurred();
            if (e) {
                jl_function_t *sprint_func = jl_get_function(jl_base_module, "sprint");
                jl_function_t *showerror_func = jl_get_function(jl_base_module, "showerror");
                jl_value_t *err_str = jl_call2(sprint_func, showerror_func, e);
                if (err_str && !jl_exception_occurred()) {
                    std::cerr << "Julia Exception: " << jl_string_ptr(err_str) << std::endl;
                }
            }
        }
        std::cerr << "Julia function returned nullptr." << std::endl;
        throw std::runtime_error("Julia function returned nullptr.");
    }
}



void LaplacianSolver::updateEdgeWeights(std::pair<int, int> edge, double weight) {
    // update the CSC graph storage
    w_edges2weights[edge] = weight;
    // Julia uses 1-based indices → stored arrays already shifted
    int u1 = edge.first + 1;
    int v1 = edge.second + 1;

    bool updated_uv = false;
    bool updated_vu = false;

    // Update edge u → v
    for (int k = colptr[v1 - 1] - 1; k < colptr[v1] - 1; ++k) {
        if (rowval[k] == u1) {
            //std::cout << "Updating edge (" << edge.first << ", " << edge.second << ") with weight " << weight << std::endl;
            nzval[k] = weight;
            updated_uv = true;
            break;
        }
    }


    // Update edge v → u (if undirected graph)
    for (int k = colptr[u1 - 1] - 1; k < colptr[u1] - 1; ++k) {
        if (rowval[k] == v1) {
            nzval[k] = weight;
            updated_vu = true;
            break;
        }
    }

    if(nzval_julia) {
        // reset the julia arrays to reflect the changes
        memcpy(jl_array_data(nzval_julia), nzval.data(), nzval.size() * sizeof(double));
    }


    if (!updated_uv || !updated_vu) {
        std::cerr << "Warning: Edge (" << edge.first << ", " << edge.second << ") not found in CSC matrix." << std::endl;
    }
}

void LaplacianSolver::NonSyncUpdateEdgeWeights(std::pair<int, int> edge, double weight) {
    // update the CSC graph storage
    w_edges2weights[edge] = weight;
    // Julia uses 1-based indices → stored arrays already shifted
    int u1 = edge.first + 1;
    int v1 = edge.second + 1;

    bool updated_uv = false;
    bool updated_vu = false;

    // Update edge u → v
    for (int k = colptr[v1 - 1] - 1; k < colptr[v1] - 1; ++k) {
        if (rowval[k] == u1) {
            //std::cout << "Updating edge (" << edge.first << ", " << edge.second << ") with weight " << weight << std::endl;
            nzval[k] = weight;
            updated_uv = true;
            break;
        }
    }


    // Update edge v → u (if undirected graph)
    for (int k = colptr[u1 - 1] - 1; k < colptr[u1] - 1; ++k) {
        if (rowval[k] == v1) {
            nzval[k] = weight;
            updated_vu = true;
            break;
        }
    }

    if (!updated_uv || !updated_vu) {
        std::cerr << "Warning: Edge (" << edge.first << ", " << edge.second << ") not found in CSC matrix." << std::endl;
    }
}

void LaplacianSolver::SyncEdgeWeights() {
    if(nzval_julia) {
        // reset the julia arrays to reflect the changes
        memcpy(jl_array_data(nzval_julia), nzval.data(), nzval.size() * sizeof(double));
    }
}