module TestZeros

# NOTE: these tests were taken from the scipy.optimize test suite

using FactCheck

include("zeros.jl")

using Zeros

#=
- f1 is a simple quadratic with roots 0 and 1
- f2 is a symmetric parabola, x^2 - 1
- f3 is a quartic polynomial with large hump in interval
- f4 is step function with a discontinuity at 1
- f5 is a hyperbola with vertical asymptote at 1
- f6 has random values positive to left of 1, negative to right

Of course these are not real problems. They just test how the 'good' solvers
behave in bad circumstances where bisection is really the best. A good solver
should not be much worse than bisection in such circumstance, while being faster
for smooth monotone sorts of functions.
=#

f1(x) = x*(x-1)
f2(x) = x^2 - 1
f3(x) = x * (x-1) * (x-2) * (x-3.)

f4(x) = x > 1 ? 1.0 + .1*x :
        x < 1 ? -1.0 + .1*x :
        0.0

f5(x) = x != 1 ? 1.0/(1. - x) : 0.0

f6(x) = x > 1 ? rand() :
        x < 1 ? -rand() :
        0.0

funcs = [f1, f2, f3, f4, f5, f6]

facts("Testing univariate root finders") do

    a = 0.5
    b = sqrt(3)

    context("testing bisect") do
        for f in funcs
            root = bisect(f, a, b; xtol=0.1e-12)
            @fact root => roughly(1.0)
        end
    end

    context("testing brent") do
        for f in funcs
            root = brent(f, a, b; xtol=0.1e-12)
            @fact root => roughly(1.0)
        end
    end

    context("testing brenth") do
        for f in funcs
            root = brenth(f, a, b; xtol=0.1e-12)
            @fact root => roughly(1.0)
        end
    end

    context("testing ridder") do
        for f in funcs
            root = ridder(f, a, b; xtol=0.1e-12)
            @fact root => roughly(1.0)
        end
    end
end

end  # module
