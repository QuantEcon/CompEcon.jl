# Tests linear basis type and constructor
module TestLin

using CompEcon
using Base.Test
using FactCheck

facts("Test Linear Basis") do

    bas1 = CompEcon.Basis(CompEcon.Lin(), [0.0,4.0], 5)

    context("test type") do
        bas2 = CompEcon.Basis(CompEcon.Lin(), [0.0,2.0,4.0], 5)
        @fact bas1.params[1].breaks == [0.0,1.0,2.0,3.0,4.0] --> true
        @fact bas2.params[1].evennum == 3 --> true
        @fact_throws ErrorException CompEcon.Basis(CompEcon.Lin(), [0.0], 5)
        @fact_throws ErrorException CompEcon.Basis(CompEcon.Lin(), [0.0,1.5,4.0], 5)
    end

    context("test constructor") do

       valbas_lin = @inferred CompEcon.evalbase(bas1.params[1])
       valbas_spl = @inferred CompEcon.evalbase(CompEcon.SplineParams([0.0,1.0,2.0,3.0,4.0],0,1))

       @fact valbas_lin == valbas_spl --> true

    end

    context("test derivative") do

    	#test derivative of basis indirectly, by fitting exponential function
		bas2 = Basis(LinParams([0,1.0],1000000))
		coeffs = funfitf(bas2,exp)
		points = rand(10000)

		d1 = funeval(coeffs,bas2,points,1)

        @fact d1 --> roughly(exp(points),atol=1e-7)

    end

end

end
