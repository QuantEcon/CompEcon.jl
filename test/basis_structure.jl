# Tests conversion among basis-structure representations
module TestBasis_Structure

using CompEcon
using Base.Test
using FactCheck


facts("Test Basis Structure Representations") do
    mb = Basis(SplineParams(linspace(-3, 3, 9), 0, 1),
               ChebParams(7, 1e-4, 10.0))

    # construct evaluation points
    X, x12 = nodes(mb)

    # construct expanded, direct, and tensor basis-structure representation
    Φ_expanded = CompEcon.BasisStructure(mb, CompEcon.Expanded(), X)
    Φ_direct = CompEcon.BasisStructure(mb, CompEcon.Direct(), X)
    Φ_tensor = CompEcon.BasisStructure(mb, CompEcon.Tensor(), x12)

    context("test convert methods") do
        @fact Φ_direct == convert(Direct, Φ_tensor) --> true
        @fact Φ_expanded == convert(Expanded, Φ_direct) --> true
        @fact Φ_expanded == convert(Expanded, Φ_tensor) --> true
        @fact ==(Φ_expanded,
                 convert(Expanded, convert(Direct, Φ_tensor))) --> true
    end

    context("test _vals_type") do
        for (TF, TM) in [(Spline, Base.SparseMatrix.SparseMatrixCSC{Float64,Int}),
                         (Lin, Base.SparseMatrix.SparseMatrixCSC{Float64,Int}),
                         (Cheb, Matrix{Float64})]
            @fact CompEcon._vals_type(TF) --> TM
            @fact CompEcon._vals_type(TF()) --> TM
        end

        @fact CompEcon._vals_type(CompEcon.BasisFamily) --> AbstractMatrix{Float64}

    end

    context("constructors") do
        nothing
    end

end  # facts


end  # module
