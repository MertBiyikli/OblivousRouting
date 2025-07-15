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



void LaplacianSolver::init(const RaeckeGraph &_g, const std::unordered_map<std::pair<int, int>, double> &_w_edges2weights, bool debug) {
    if(initted)
        return;

    m_graph = _g;
    for(const auto& [edge, weight] : _w_edges2weights) {
        w_edges2weights[edge] = weight;
    }
    m = m_graph.getNumEdges();
    n = m_graph.getNumNodes();

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
    bvec[0]=1;
    bvec[n-1]=-1;

    std::cout << "bvec to send to Julia:\n" << bvec.transpose() << std::endl;


    buildCSC(_w_edges2weights);



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


Eigen::VectorXd LaplacianSolver::solve(const Eigen::VectorXd &b, bool debug) {

    if(!initted) {
        throw std::runtime_error("LaplacianSolver not initialized. Call init() first.");
    }
    if(b.size() != n) {
        throw std::runtime_error("Input vector size does not match the number of nodes in the graph.");
    }

    // init the bvec if the dimensions are matching
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


    jl_value_t* args[6] = {
            arg1,
            arg2,
            (jl_value_t *) colptr_julia,
            (jl_value_t *) rowval_julia,
            (jl_value_t *) nzval_julia,
            (jl_value_t *) bvec_julia
    };

    jl_array_t* result = (jl_array_t *)jl_call(func, args, 6);

    if(debug) {
        if (jl_exception_occurred()) {
            jl_value_t *e = jl_exception_occurred();
            jl_function_t *sprint_func = jl_get_function(jl_base_module, "sprint");
            jl_function_t *showerror_func = jl_get_function(jl_base_module, "showerror");
            jl_value_t *err_str = jl_call2(sprint_func, showerror_func, e);
            if (err_str && !jl_exception_occurred()) {
                std::cerr << "Julia Exception: " << jl_string_ptr(err_str) << std::endl;
            }
            throw std::runtime_error("Julia call failed.");
        }
    }
    // Convert the result to Eigen::VectorXd
    double* data_ptr = (double*)jl_array_data(result);
    size_t len = jl_array_len(result);

    Eigen::VectorXd x(len);
    for (size_t i = 0; i < len; ++i) {
        x[i] = data_ptr[i];
    }

    return x;

}



void LaplacianSolver::buildCSC(const std::unordered_map<std::pair<int, int>, double> &edges) {
    // Prepare edge list
    std::vector<std::tuple<int, int, double>> edgeList;
    edgeList.reserve(edges.size() * 2);

    for (const auto &[edge, weight] : edges) {
        if (edge.first == edge.second) continue; // skip self-loops

        // add both directions to ensure symmetry
        edgeList.emplace_back(edge.first, edge.second, weight);
        edgeList.emplace_back(edge.second, edge.first, weight);
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

    // Debug prints
    std::cout << "colptr:";
    for (auto v : colptr) std::cout << " " << v;
    std::cout << std::endl;

    std::cout << "rowval:";
    for (auto v : rowval) std::cout << " " << v;
    std::cout << std::endl;

    std::cout << "nzval:";
    for (auto v : nzval) std::cout << " " << v;
    std::cout << std::endl;
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
            std::cout << "Updating edge (" << edge.first << ", " << edge.second << ") with weight " << weight << std::endl;
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
    std::cout << "Updated nzval:";
    for (auto v : nzval) std::cout << " " << v;
    std::cout << std::endl;

    if (!updated_uv || !updated_vu) {
        std::cerr << "Warning: Edge (" << edge.first << ", " << edge.second << ") not found in CSC matrix." << std::endl;
    }
}