# Tests constructors and subtype relationships
module TestTypeSystem

using CompEcon
using Base.Test
using FactCheck

sp = SplineParams(linspace(0, 5, 10), 0, 3)
cp = ChebParams(6, 2.0, 5.0)
lp = LinParams(linspace(0, 5, 10), 0)

facts("Test subtype structure") do
    context("Basis families") do
        for T in [Cheb, Lin, Spline]
            @fact T <: BasisFamily --> true
        end
    end

    context("params") do
        for T in [ChebParams, SplineParams, LinParams]
            @fact T <: BasisParams --> true
        end
    end

    context("basis structure representations") do
        for T in [Tensor, Direct, Expanded]
            @fact T <: AbstractBasisStructureRep --> true
        end
    end
end

facts("test _param method") do
    for (T, TP) in [(Cheb, ChebParams), (Lin, LinParams), (Spline, SplineParams)]
        @fact CompEcon._param(T) --> TP
        @fact CompEcon._param(T()) --> TP
    end
end

facts("Test original API compat") do

    s = Spline()
    c = Cheb()
    l = Lin()

    context("old_name for basis family") do
        @fact old_name(s) --> :spli
        @fact old_name(c) --> :cheb
        @fact old_name(l) --> :lin
    end

    context("old_name for params") do
        @fact old_name(sp) --> :spli
        @fact old_name(cp) --> :cheb
        @fact old_name(lp) --> :lin
    end

    context("old_params") do
        @fact CompEcon.old_params(sp) --> Any[sp.breaks, sp.evennum, sp.k]
        @fact CompEcon.old_params(cp) --> Any[cp.n, cp.a, cp.b]
        @fact CompEcon.old_params(lp) --> Any[lp.breaks, lp.evennum]
    end

    context("convert(Basis, Dict) and revert(Basis) --> Dict") do
        b = Basis(sp, cp)
        d = revert(b)
        @fact convert(Basis, d) --> b
    end

end

facts("test more Param methods") do

    context("Test extra outer constructors") do
        @fact sp --> SplineParams(10, 0, 5, 3)
        @fact lp --> LinParams(10, 0, 5)
    end

    context("test equality of params") do
        @fact ==(sp, lp) --> false
        @fact ==(sp, cp) --> false
        @fact ==(cp, lp) --> false

        @fact ==(sp, sp) --> true
        @fact ==(cp, cp) --> true
        @fact ==(lp, lp) --> true
    end
end

end  # module TestTypeSystem
