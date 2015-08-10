# Tests conversion among basis-structure representations
module TestInterp

# include instead of `using` so it updates when I reload the module
include(Pkg.dir("CompEcon", "src", "CompEcon.jl"))
using Base.Test
using FactCheck

facts("check interp.jl") do

    # construct basis
    basis = CompEcon.Basis(CompEcon.SplineParams(15, -1, 1, 1),CompEcon.SplineParams(20, -5, 2, 3))

    # get nodes
    X, x12 = CompEcon.nodes(basis)

    # function to interpolate
    f(x1, x2) = cos(x1) ./ exp(x2)
    f(X::Matrix) = f(X[:, 1], X[:, 2])

    # function at nodes
    y = f(X)

    # benchmark coefficients
    c_direct, bs_direct = CompEcon.funfitxy(basis, X, y)

    context("test funfitxy for tensor and direct agree on coefs") do
        c_tensor, bs_tensor = CompEcon.funfitxy(basis, x12, y)
        @fact c_tensor --> roughly(c_direct; atol=1e-13)
    end

    context("test funfitf") do
        c = CompEcon.funfitf(basis,f)

        @fact c --> roughly(c_direct; atol=1e-15)
    end

    context("test funeval methods") do
        #single point
        sp = CompEcon.funeval(c_direct,basis,X[5,:])[1]
        @fact sp --> roughly(y[5]; atol=1e-15)

        #multiple points using tensor directly
        mp = CompEcon.funeval(c_direct,basis,x12)
        @fact mp --> roughly(y; atol=1e-15)

        #multiple points using direct (which in turn uses tensor)
        mp = CompEcon.funeval(c_direct,basis,X)
        @fact mp --> roughly(y; atol=1e-15)

        #multiple points giving basis in direct form
        Phidirect = CompEcon.BasisStructure(basis,CompEcon.Direct(),X)
        mpd = CompEcon.funeval(c_direct,Phidirect)
        @fact mpd --> roughly(y; atol=1e-15)

        #multiple points giving basis in expanded form
        Phiexp = Base.convert(CompEcon.Expanded,Phidirect)
        mpe = CompEcon.funeval(c_direct,Phiexp)
        @fact mpe --> roughly(y; atol=1e-15)

    end

    context("test interpoland methods") do
        # (Basis,BasisStructure,..)
        Phidirect = CompEcon.BasisStructure(basis,CompEcon.Direct(),X)
        intp1 = CompEcon.Interpoland(basis,Phidirect,y)
        @fact CompEcon.evaluate(intp1,X) --> roughly(y; atol=1e-15)
        # (Basis,Array,..)
        intp2 = CompEcon.Interpoland(basis,X,y)
        @fact CompEcon.evaluate(intp2,X) --> roughly(y; atol=1e-15)
        # (BasisParams,Function)
        intp3 = CompEcon.Interpoland(basis,f)
        @fact CompEcon.evaluate(intp3,X) --> roughly(y; atol=1e-15)
    end

    # TODO: call show on an interpoland instance to get coverage for writemime

end  # facts

end  # module
