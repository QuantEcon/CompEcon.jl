# Tests conversion among basis-structure representations
module TestInterp

# include instead of `using` so it updates when I reload the module
include(Pkg.dir("CompEcon", "src", "CompEcon.jl"))
using Base.Test
using FactCheck

facts("check interp.jl") do

    # construct basis
    basis = CompEcon.Basis(CompEcon.SplineParams(15, -1, 1, 1),
                           CompEcon.SplineParams(20, -5, 2, 3))

    # get nodes
    X, x12 = CompEcon.nodes(basis)

    # function to interpolate
    f(x1, x2) = cos(x1) ./ exp(x2)
    f(X::Matrix) = f(X[:, 1], X[:, 2])

    # function at nodes
    y = f(X)

    context("test funfitxy for tensor and direct agree on coefs") do
        c_tensor, bs_tensor = CompEcon.funfitxy(basis, x12, y)
        c_direct, bs_direct = CompEcon.funfitxy(basis, X, y)

        @fact c_tensor --> rougly(c_direct; atol=1e-13)
    end

    # TODO: call show on an interpoland instance to get coverage for writemime

end  # facts

end  # module
