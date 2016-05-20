# Tests spline basis evaluation function
module TestSpline

using CompEcon
using Base.Test
using FactCheck

facts("Test Spline Basis Evaluation") do

    # constuct specific case that we can compute by hand

    n = 5
    a = 0
    b = 1
    params = SplineParams(linspace(a,b,n),0,1)

    x = rand(10000)

    basestr = @inferred CompEcon.evalbase(params, x)
    base = full(basestr)

    # vector of basis functions

    function linBspline(u)

        if u >= 0 && u < 0.25

            return [1-4*u,4*u,0,0,0]

        elseif u >= 0.25 && u < 0.5

            return [0,2*(1 - 2*u),4*u-1,0,0]

        elseif u >= 0.5 && u < 0.75

            return [0,0,3 - 4*u,2*(2*u - 1),0]

        else

            return [0,0,0,4*(1 - u),4*u-3]

        end

    end

    manualbase = Array(Float64,length(x),n)

    for i in 1:length(x)

        manualbase[i,:] = linBspline(x[i])

    end

    context("test evalbase with linear B spline") do
        @fact manualbase == base --> true

    end

end

end
