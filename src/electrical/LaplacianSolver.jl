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
    bvec::Vector{Float64}
)

    @show typeof(n)
    @show typeof(m)
    @show typeof(colptr)
    @show typeof(rowval)
    @show typeof(nzval)
    @show typeof(bvec)

    println("length(bvec) = ", length(bvec))
    println("eltype(bvec) = ", eltype(bvec))


    println("sum(bvec) = ", sum(bvec))

    bvec = bvec .- mean(bvec)
    println("sum(bvec) = ", sum(bvec))
#=

    U = zeros(Float64, n)
    L = SparseMatrixCSC(n, m, colptr, rowval, nzval)
    solver = approxchol_sddm(L, verbose=false)
    U = solver(Vector(bvec), verbose=false)
=#
    # Example where we output the adjacaency matrix
    a = ring_graph(50)
    @show typeof(a)
    @show size(a)
    @show eltype(a)
    @show a[1:5, 1:5]

    # Create a sparse matrix from the provided data
    A = SparseMatrixCSC(n, m, colptr, rowval, nzval)
    @show typeof(A)
    @show size(A)
    @show eltype(A)
    @show A[1:5, 1:5]
    sol = approxchol_lap(A)
    x = sol(bvec, tol=1e-1, verbose=true)
    return x
end

end # module
