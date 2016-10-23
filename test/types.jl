sp = SplineParams(linspace(0, 5, 10), 0, 3)
cp = ChebParams(6, 2.0, 5.0)
lp = LinParams(linspace(0, 5, 10), 0)

@testset "Test subtype structure" begin
    @testset "Basis families" begin
        for T in [Cheb, Lin, Spline]
            @test T <: BasisFamily
        end
    end

    @testset "params" begin
        for T in [ChebParams, SplineParams, LinParams]
            @test T <: BasisParams
        end
    end

    @testset "basis structure representations" begin
        for T in [Tensor, Direct, Expanded]
            @test T <: AbstractBasisStructureRep
        end
    end
end

@testset "test _param method" begin
    for (T, TP) in [(Cheb, ChebParams), (Lin, LinParams), (Spline, SplineParams)]
        @test CompEcon._param(T) == TP
        @test CompEcon._param(T()) == TP
    end
end

@testset "Test original API compat" begin

    s = Spline()
    c = Cheb()
    l = Lin()

    @testset "old_name for basis family" begin
        @test Original.old_name(s) == :spli
        @test Original.old_name(c) == :cheb
        @test Original.old_name(l) == :lin
    end

    @testset "old_name for params" begin
        @test Original.old_name(sp) == :spli
        @test Original.old_name(cp) == :cheb
        @test Original.old_name(lp) == :lin
    end

    @testset "old_params" begin
        @test Original.old_params(sp) == Any[sp.breaks, sp.evennum, sp.k]
        @test Original.old_params(cp) == Any[cp.n, cp.a, cp.b]
        @test Original.old_params(lp) == Any[lp.breaks, lp.evennum]
    end

    @testset "convert(Basis, Dict) and revert(Basis) --> Dict" begin
        b = Basis(sp, cp)
        d = Original.revert(b)
        @test convert(Basis, d) == b
    end

end

@testset "test more Param methods" begin

    @testset "Test extra outer constructors" begin
        @test sp == SplineParams(10, 0, 5, 3)
        @test lp == LinParams(10, 0, 5)
    end

    @testset "test equality of params" begin
        @test (  ==(sp, lp) ) == false
        @test (  ==(sp, cp) ) == false
        @test (  ==(cp, lp) ) == false

        @test   ==(sp, sp)
        @test   ==(cp, cp)
        @test   ==(lp, lp)
    end
end
