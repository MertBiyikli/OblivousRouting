module LaplacianSolver

import Pkg
Pkg.activate(joinpath(@__DIR__, "../../../"))

@info "Activated Julia project:" Base.active_project()

using SparseArrays
using Laplacians
using Statistics

export laplacian_solve

function laplacian_solve(
    n::Int,
    m::Int,
    colptr::Vector{Int32},
    rowval::Vector{Int32},
    nzval::Vector{Float64},
    bvec::Vector{Float64},
    debug::Bool = false,
)

    if debug
        @show typeof(n)
        @show typeof(m)
        @show typeof(colptr)
        @show typeof(rowval)
        @show typeof(nzval)
        @show typeof(bvec)

        println("JULIA: colptr length = ", length(colptr))
        println("JULIA: rowval length = ", length(rowval))
        println("JULIA: nzval length = ", length(nzval))
        println("JULIA: first few rowval = ", rowval[1:min(end, 10)])
        println("JULIA: first few nzval = ", nzval[1:min(end, 10)])

        println("length(bvec) = ", length(bvec))
        println("eltype(bvec) = ", eltype(bvec))
        println("sum(bvec) = ", sum(bvec))

        @assert length(rowval) == length(nzval) "rowval and nzval must be the same length"
        @assert colptr[end] - 1 == length(rowval) "colptr mismatch: total nnz ≠ rowval/nzval length"
    end

    bvec = bvec .- mean(bvec)

    # Create a sparse matrix from the provided data
    A = SparseMatrixCSC(n, n, colptr, rowval, nzval)
    if debug
        @show size(A)
        @show nnz(A)
        @show typeof(A)
    end

    sol = approxchol_lap(A)
    x = sol(bvec, tol=1e-10, verbose=false)
    return x
end

end # module
