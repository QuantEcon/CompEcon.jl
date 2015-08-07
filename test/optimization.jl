# Tests golden method
module TestOptimization

using CompEcon
using Base.Test
using FactCheck


facts("Test Optimization") do

	# test 1D method

	f(x) = -(x^2 + x) # the function has a global maximum at x = -.5
	a = -3
	b = 5

	context("Golden Method, 1D") do
        xstar,fstar = golden_method(f,a,b)
        @fact -.5 --> roughly(xstar, atol = 1e-15)
    	@fact .25 --> roughly(fstar, atol = 1e-15)
    end

	# test multiD method

	g(x) = [f(x[1]),-abs(x[2])] # the 1st function has a global maximum at x = -.5, the 2nd has a global maximum at x = 0
	a = [-3,-5]
	b = [5,2]

	context("Golden Method, multiD") do
        xstar,fstar = golden_method(g,a,b)
        @fact [-.5,0] --> roughly(xstar, atol = 1e-15)
        @fact [.25,0] --> roughly(fstar, atol = 1e-15)
    end

end

end
