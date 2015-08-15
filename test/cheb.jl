# Tests Cheb basis
module TestCheb

# include instead of `using` so it updates when I reload the module
include(Pkg.dir("CompEcon", "src", "CompEcon.jl"))
using Base.Test
using FactCheck

params = CompEcon.ChebParams(7, 1e-4, 10.0)

nod = CompEcon.nodes(params)
nod_end = CompEcon.nodes(params,1) # extended nodes to endpoints
nod_l = CompEcon.nodes(params,3) # Lobatto nodes

facts("test Cheb") do

    context("test nodes") do
        @fact  (nod[1] == params.a[1],nod[end] == params.b[1]) --> (false,false)
        @fact  [nod_end[1],nod_end[end]] --> roughly([params.a[1],params.b[1]], atol = 1e-15)
        @fact  [nod_l[1],nod_l[end]] --> roughly([params.a[1],params.b[1]], atol = 1e-15)
    end

    context("test derivative/integral") do
        deriv, ord_d, a_d, b_d, holder_d = CompEcon.derivative_op(params,1)
        int, ord_i, a_i, b_i, holder_i = CompEcon.derivative_op(params,-1)

        @fact  (ord_d == params.n[1]-1,ord_i == params.n[1]+1) --> (true,true)
        @fact  [a_d,b_d] --> roughly([params.a[1],params.b[1]], atol = 1e-15)
        @fact  [a_i,b_i] --> roughly([params.a[1],params.b[1]], atol = 1e-15)
    end

    context("test evalbase") do
        B, x = CompEcon.evalbase(params,nod,0,0)
        B_end, x_end = CompEcon.evalbase(params,nod_end,0,1) # extended nodes to endpoints
    end

    # TODO: call show to get coverage for writemime

end  # facts

end  # module
