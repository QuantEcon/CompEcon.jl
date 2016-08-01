# Tests Cheb basis
module TestCheb

# include instead of `using` so it updates when I reload the module
include(joinpath(dirname(@__FILE__), "..", "src", "CompEcon.jl"))
using Base.Test
using FactCheck

params = CompEcon.ChebParams(7, 1e-4, 10.0)

nod = CompEcon.nodes(params)
nod_end = CompEcon.nodes(params,1) # extended nodes to endpoints
nod_l = CompEcon.nodes(params,3) # Lobatto nodes

function manualeval(x,a,b,n)

    z = 2*(x-a)/(b-a)-1
    m = length(x)

    bas = Array(Float64,m,n)
    bas[:,1:2] = [ones(m,1) z]

    for i in 3:n

        bas[:,i] = 2*z.*bas[:,i-1]-bas[:,i-2]

    end

    return bas

end

mB = manualeval(nod,params.a[1],params.b[1],params.n[1])

#=
function manualderiv1(a,b,n)

    x = CompEcon.nodes(CompEcon.ChebParams(n-1,a,b))
    z = (2*(x-a)/(b-a)-1)''

    D = [zeros(length(z),1) ones(length(z),1) 4*z 12*z.^2-3 32*z.^3-16*z 80*z.^4-60*z.^2+5 192*z.^5-192*z.^3+36*z]

    return D

end

function manualint1(a,b,n)

    x = CompEcon.nodes(CompEcon.ChebParams(n+1,a,b))
    z = (2*(x-a)/(b-a)-1)''

    I = [z .5*z.^2 2/3*z.^3-z z.^4-3/2*z.^2 8/5*z.^5-8/3*z.^3+z 16/6*z.^6-5*z.^4+2.5*z.^2 32/7*z.^7-48/5*z.^5+6*z.^3-z]

    return I

end

D = manualderiv1(nod,params.a[1],params.b[1])
I = manualint1(nod,params.a[1],params.b[1])
=#

facts("test Cheb") do

    context("test nodes") do
        @fact  (nod[1] == params.a[1],nod[end] == params.b[1]) --> (false,false)
        @fact  [nod_end[1],nod_end[end]] --> roughly([params.a[1],params.b[1]], atol = 1e-15)
        @fact  [nod_l[1],nod_l[end]] --> roughly([params.a[1],params.b[1]], atol = 1e-15)
    end

    context("test derivative/integral") do
        deriv, params_d = CompEcon.derivative_op(params,1)
        int, params_i = CompEcon.derivative_op(params,-1)

        @fact  params_d.n == params.n[1]-1 --> true
        @fact  params_i.n == params.n[1]+1 --> true
        @fact  [params_d.a,params_d.b] --> roughly([params.a[1],params.b[1]], atol = 1e-15)
        @fact  [params_i.a,params_i.b] --> roughly([params.a[1],params.b[1]], atol = 1e-15)
    end

    context("test evalbase") do
        B = @inferred CompEcon.evalbase(params,nod,0,0)
        @fact  B --> roughly(mB, atol = 1e-14)
    end

end  # facts

end  # module
