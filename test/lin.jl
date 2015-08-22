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

       valbas_lin = CompEcon.evalbase(bas1.params[1])
       valbas_spl = CompEcon.evalbase(CompEcon.SplineParams([0.0,1.0,2.0,3.0,4.0],0,1))

       @fact valbas_lin == valbas_spl --> true

    end

end

end
