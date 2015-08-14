# Tests Cheb basis
module TestCheb

# include instead of `using` so it updates when I reload the module
include(Pkg.dir("CompEcon", "src", "CompEcon.jl"))
using Base.Test
using FactCheck

params = CompEcon.ChebParams(7, 1e-4, 10.0)

basis = CompEcon.Basis(params)

nodes = nodes(params)
nodes_end = nodes(params,1) # extended nodes to endpoints
nodes_l = nodes(params,3) # Lobatto nodes

facts("test Cheb") do

    context("test derivative/integral") do
        
        deriv, ord_d, a_d, b_d, holder_d = derivative_op(params)
        int, ord_i, a_i, b_i, holder_i = derivative_op(params,-1)

        @fact  --> 
    end


    # TODO: call show to get coverage for writemime

end  # facts

end  # module
