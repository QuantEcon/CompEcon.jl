# Tests conversion among basis-structure representations

@testset "Test Basis Structure Representations" begin
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

    @testset "test standard Base methods" begin
        # test == and ndims, multiD
        for Φ in (Φ_expanded, Φ_direct, Φ_tensor)
            @test  Φ  ==  Φ
            @test ndims(Φ) == 3
        end

        for Φ in (Φ_direct, Φ_tensor)
            @test ( Φ_expanded  ==  Φ ) == false
        end

        # test == and ndims, 1D
        for Φ in (Φ_expanded_1d, Φ_direct_1d)
            @test  Φ  ==  Φ
            @test ndims(Φ) == 1
        end

        @test  Φ_expanded_1d  ==  Φ_direct_1d

    end

    @testset "test convert methods" begin
        @test  Φ_direct  ==  convert(Direct, Φ_tensor)
        @test  Φ_expanded  ==  convert(Expanded, Φ_direct)
        @test  Φ_expanded  ==  convert(Expanded, Φ_tensor)
        @test ==(Φ_expanded,
                 convert(Expanded, convert(Direct, Φ_tensor)))


        # test the no-op covert method
        @test convert(Expanded, Φ_expanded) == Φ_expanded
        @test convert(Direct, Φ_direct) == Φ_direct
        @test convert(Tensor, Φ_tensor) == Φ_tensor

        # test that we can't do expanded -> (direct|tensor) or direct -> tensor
        @test_throws MethodError convert(Direct, Φ_expanded)
        @test_throws MethodError convert(Tensor, Φ_expanded)
        @test_throws MethodError convert(Tensor, Φ_direct)
    end

    @testset "test internal tools" begin
        ## test _vals_type
        for (TF, TM) in [(Spline, SparseMatrixCSC{Float64,Int}),
                         (Lin, SparseMatrixCSC{Float64,Int}),
                         (Cheb, Matrix{Float64})]
            @test CompEcon._vals_type(TF) == TM
            @test CompEcon._vals_type(TF()) == TM
        end

        @test CompEcon._vals_type(CompEcon.BasisFamily) == AbstractMatrix{Float64}

        ## test _checkx
        # create test data
        xm = rand(1, 2)
        xv = rand(2)
        xvv = [rand(2) for i=1:2]

        @test CompEcon._checkx(2, xm) == xm
        @test CompEcon._checkx(2, xv) == reshape(xv, 1, 2)
        @test CompEcon._checkx(2, xvv) == xvv
        @test CompEcon._checkx(1, xv) == xv

        @test_throws ErrorException CompEcon._checkx(2, rand(1, 3))
        @test_throws ErrorException CompEcon._checkx(2, rand(3))
        @test_throws ErrorException CompEcon._checkx(2, [rand(2) for i=1:3])

        ## test check_convert
        for Φ in (Φ_expanded, Φ_direct, Φ_tensor)
            @test CompEcon.check_convert(Φ, zeros(1, 3)) == (3, 1, 3)
            @test_throws ErrorException CompEcon.check_convert(Φ, zeros(1, 2))
            @test_throws ErrorException CompEcon.check_convert(Φ, fill(-1, 1, 3))
        end

        ## test check_basis_structure
        @testset "test check_basis_structure" begin
            order1 = zeros(Int, 1, 2)
            out1 = CompEcon.check_basis_structure(2, xm, order1)
            @test out1 == (1, order1, order1, [1 1], xm)

            # test when order::Int
            @test CompEcon.check_basis_structure(2, xm, 0) == out1

            # test the reshape(order, 1, N) branch (isa(order, Vector))
            @test CompEcon.check_basis_structure(2, xm, [0, 0]) == out1

            # check N=1 --> order = fill(order, 1, 1) branch
            out_1d = CompEcon.check_basis_structure(1, xv, 0)
            order_1d_out = fill(0, 1, 1)
            @test out_1d == (1, order_1d_out, order_1d_out, fill(1, 1, 1), xv)

            # check m > 1 branch
            order_m2 = [0 0; 1 1; 0 -1]
            out_m2 = CompEcon.check_basis_structure(2, xm, order_m2)
            @test out_m2 == (3, order_m2, [0 -1], [2 3], xm)


            order2 = zeros(Int, 1, 5)  # should throw error
            @test_throws ErrorException CompEcon.check_basis_structure(2, X,
                                                                       order2)
        end


    end

    @testset "constructors" begin
        nothing
    end

    @testset "Test from PR #25" begin
        basisμ = CompEcon.Basis(Cheb, 20, 0.0, 1.0)
        basisσ = CompEcon.Basis(Cheb, 20, 0.0, 1.0)
        basis = CompEcon.Basis(basisμ, basisσ)
        S, (μs, σs) = CompEcon.nodes(basis)
        bs = CompEcon.BasisStructure(basis, CompEcon.Expanded(), S, [0 2])
        @test isa(bs, CompEcon.BasisStructure{Expanded}) == true
    end


    # call show (which calls writemime) so we can get 100% coverage
    @testset "Printing" begin
        iob = IOBuffer()
        show(iob, Φ_tensor)
    end

end  # testset
