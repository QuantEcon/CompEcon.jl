cases = ["B-splines","Cheb","Lin"]

# construct basis
holder = (CompEcon.Basis(CompEcon.SplineParams(15,-1,1,1),CompEcon.SplineParams(20,-5,2,3)),
    CompEcon.Basis(CompEcon.ChebParams(15,-1,1),CompEcon.ChebParams(20,-5,2)),
    CompEcon.Basis(CompEcon.LinParams(15,-1,1),CompEcon.LinParams(20,-5,2)))

@testset for (i, case) in enumerate(cases)
    basis = holder[i]

    # get nodes
    X, x12 = CompEcon.nodes(basis)

    # function to interpolate
    f(x1, x2) = cos(x1) ./ exp(x2)
    f(X::Matrix) = f(X[:, 1], X[:, 2])

    # function at nodes
    y = f(X)

    # benchmark coefficients
    c_direct, bs_direct = CompEcon.funfitxy(basis, X, y)

    @testset "test funfitxy for tensor and direct agree on coefs" begin
        c_tensor, bs_tensor = CompEcon.funfitxy(basis, x12, y)
        @test maxabs(c_tensor -  c_direct) <=  1e-12
    end

    @testset "test funfitf" begin
        c = CompEcon.funfitf(basis,f)
        @test maxabs(c -  c_direct) <=  1e-12
    end

    @testset "test funeval methods" begin
        # single point
        sp = CompEcon.funeval(c_direct,basis,X[5:5,:])[1]
        @test maxabs(sp -  y[5]) <= 1e-12

        # multiple points using tensor directly
        mp = CompEcon.funeval(c_direct,basis,x12)
        @test maxabs(mp -  y) <=  1e-12

        # multiple points using direct
        mp = CompEcon.funeval(c_direct,basis,X)
        @test maxabs(mp -  y) <=  1e-12

        # multiple points giving basis in direct form
        mpd = CompEcon.funeval(c_direct,bs_direct)
        @test maxabs(mpd -  y) <=  1e-12

        # multiple points giving basis in expanded form
        Phiexp = Base.convert(CompEcon.Expanded,bs_direct)
        mpe = CompEcon.funeval(c_direct,Phiexp)
        @test maxabs(mpe -  y) <=  1e-12

    end

    @testset "test interpoland methods" begin
        # (Basis,BasisStructure,..)
        intp1 = CompEcon.Interpoland(basis,bs_direct,y)
        @test maxabs(CompEcon.evaluate(intp1,X) - y) <= 1e-12
        # (Basis,Array,..)
        intp2 = CompEcon.Interpoland(basis,X,y)
        @test maxabs(CompEcon.evaluate(intp2,X) - y) <= 1e-12
        # (BasisParams,Function)
        intp3 = CompEcon.Interpoland(basis,f)
        @test maxabs(CompEcon.evaluate(intp3,X) - y) <= 1e-12
    end

    @testset "Printing" begin
        iob = IOBuffer()
        show(iob, CompEcon.Interpoland(basis, bs_direct, y))
    end

    # TODO: call show on an interpoland instance to get coverage for writemime

end # testset
