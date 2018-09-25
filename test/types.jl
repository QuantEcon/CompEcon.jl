sp = SplineParams(range(0, stop = 5, length = 10), 0, 3)
cp = ChebParams(6, 2.0, 5.0)
lp = LinParams(range(0, stop = 5, length = 10), 0)

@testset "Test original API compat" begin

    s = Spline()
    c = Cheb()
    l = Lin()

    @testset "old_name for basis family" begin
        @test CompEcon.old_name(s) == :spli
        @test CompEcon.old_name(c) == :cheb
        @test CompEcon.old_name(l) == :lin
    end

    @testset "old_name for params" begin
        @test CompEcon.old_name(sp) == :spli
        @test CompEcon.old_name(cp) == :cheb
        @test CompEcon.old_name(lp) == :lin
    end

    @testset "old_params" begin
        @test CompEcon.old_params(sp) == Any[sp.breaks, sp.evennum, sp.k]
        @test CompEcon.old_params(cp) == Any[cp.n, cp.a, cp.b]
        @test CompEcon.old_params(lp) == Any[lp.breaks, lp.evennum]
    end

    @testset "convert(Basis, Dict) and revert(Basis) --> Dict" begin
        b = Basis(sp, cp)
        d = CompEcon.revert(b)
        @test convert(Basis, d) == b
    end

end
