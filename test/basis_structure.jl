# Tests conversion among basis-structure representations
module TestBasisStructure

using CompEcon
using Base.Test
using FactCheck


facts("Test Basis Structure Representations") do
    mb = Basis(SplineParams(linspace(-3, 3, 9), 0, 1),
               ChebParams(7, 1e-4, 10.0),
               LinParams(linspace(-4, 2, 11), 0))

    # construct evaluation points
    X, x123 = nodes(mb)

    # construct expanded, direct, and tensor basis-structure representation
    Φ_expanded = CompEcon.BasisStructure(mb, CompEcon.Expanded(), X)
    Φ_direct = CompEcon.BasisStructure(mb, CompEcon.Direct(), X)
    Φ_tensor = CompEcon.BasisStructure(mb, CompEcon.Tensor(), x123)

    # construct expanded, direct, and tensor basis-structure representation with 1D basis
    Φ_expanded_1d = CompEcon.BasisStructure(mb[1], CompEcon.Expanded(), X[:,1])
    Φ_direct_1d = CompEcon.BasisStructure(mb[1], CompEcon.Direct(), X[:,1])

    context("test standard Base methods") do
        # test == and ndims, multiD
        for Φ in (Φ_expanded, Φ_direct, Φ_tensor)
            @fact Φ == Φ --> true
            @fact ndims(Φ) --> 3
        end

        for Φ in (Φ_direct, Φ_tensor)
            @fact Φ_expanded == Φ --> false
        end

        # test == and ndims, 1D
        for Φ in (Φ_expanded_1d, Φ_direct_1d)
            @fact Φ == Φ --> true
            @fact ndims(Φ) --> 1
        end

        @fact Φ_expanded_1d == Φ_direct_1d --> true

    end

    context("test convert methods") do
        @fact Φ_direct == convert(Direct, Φ_tensor) --> true
        @fact Φ_expanded == convert(Expanded, Φ_direct) --> true
        @fact Φ_expanded == convert(Expanded, Φ_tensor) --> true
        @fact ==(Φ_expanded,
                 convert(Expanded, convert(Direct, Φ_tensor))) --> true


        # test the no-op covert method
        @fact convert(Expanded, Φ_expanded) --> Φ_expanded
        @fact convert(Direct, Φ_direct) --> Φ_direct
        @fact convert(Tensor, Φ_tensor) --> Φ_tensor

        # test that we can't do expanded -> (direct|tensor) or direct -> tensor
        @fact_throws MethodError convert(Direct, Φ_expanded)
        @fact_throws MethodError convert(Tensor, Φ_expanded)
        @fact_throws MethodError convert(Tensor, Φ_direct)
    end

    context("test internal tools") do
        ## test _vals_type
        for (TF, TM) in [(Spline, Base.SparseMatrix.SparseMatrixCSC{Float64,Int}),
                         (Lin, Base.SparseMatrix.SparseMatrixCSC{Float64,Int}),
                         (Cheb, Matrix{Float64})]
            @fact CompEcon._vals_type(TF) --> TM
            @fact CompEcon._vals_type(TF()) --> TM
        end

        @fact CompEcon._vals_type(CompEcon.BasisFamily) --> AbstractMatrix{Float64}

        ## test _checkx
        # create test data
        xm = rand(1, 2)
        xv = rand(2)
        xvv = [rand(2) for i=1:2]

        @fact CompEcon._checkx(2, xm) --> xm
        @fact CompEcon._checkx(2, xv) --> reshape(xv, 1, 2)
        @fact CompEcon._checkx(2, xvv) --> xvv
        @fact CompEcon._checkx(1, xv) --> xv

        @fact_throws ErrorException CompEcon._checkx(2, rand(1, 3))
        @fact_throws ErrorException CompEcon._checkx(2, rand(3))
        @fact_throws ErrorException CompEcon._checkx(2, [rand(2) for i=1:3])

        ## test check_convert
        for Φ in (Φ_expanded, Φ_direct, Φ_tensor)
            @fact CompEcon.check_convert(Φ, zeros(1, 3)) --> (3, 1, 3)
            @fact_throws ErrorException CompEcon.check_convert(Φ, zeros(1, 2))
            @fact_throws ErrorException CompEcon.check_convert(Φ, fill(-1, 1, 3))
        end

        ## test check_basis_structure
        context("test check_basis_structure") do
            order1 = zeros(Int, 1, 2)
            out1 = CompEcon.check_basis_structure(2, xm, order1)
            @fact out1 --> (1, order1, order1, [1 1], xm)

            # test when order::Int
            @fact CompEcon.check_basis_structure(2, xm, 0) --> out1

            # test the reshape(order, 1, N) branch (isa(order, Vector))
            @fact CompEcon.check_basis_structure(2, xm, [0, 0]) --> out1

            # check N=1 --> order = fill(order, 1, 1) branch
            out_1d = CompEcon.check_basis_structure(1, xv, 0)
            order_1d_out = fill(0, 1, 1)
            @fact out_1d --> (1, order_1d_out, order_1d_out, fill(1, 1, 1), xv)

            # check m > 1 branch
            order_m2 = [0 0; 1 1; 0 -1]
            out_m2 = CompEcon.check_basis_structure(2, xm, order_m2)
            @fact out_m2 --> (3, order_m2, [0 -1], [2 3], xm)


            order2 = zeros(Int, 1, 5)  # should throw error
            @fact_throws ErrorException CompEcon.check_basis_structure(2, X,
                                                                       order2)
        end


    end

    context("constructors") do
        nothing
    end

    context("Test from PR #25") do
        basisμ = CompEcon.Basis(Cheb, 20, 0.0, 1.0)
        basisσ = CompEcon.Basis(Cheb, 20, 0.0, 1.0)
        basis = CompEcon.Basis(basisμ, basisσ)
        S, (μs, σs) = CompEcon.nodes(basis)
        bs = CompEcon.BasisStructure(basis, CompEcon.Expanded(), S, [0 2])
        @fact isa(bs, CompEcon.BasisStructure{Expanded}) --> true
    end


    # call show (which calls writemime) so we can get 100% coverage
    context("Printing") do
        iob = IOBuffer()
        show(iob, Φ_tensor)
    end

end  # facts


end  # module
