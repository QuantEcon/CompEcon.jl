@testset "Test Linear Basis" begin

    bas1 = CompEcon.Basis(CompEcon.Lin(), [0.0,4.0], 5)

    @testset "test type" begin
        bas2 = CompEcon.Basis(CompEcon.Lin(), [0.0,2.0,4.0], 5)
        @test  bas1.params[1].breaks  ==  [0.0,1.0,2.0,3.0,4.0]
        @test  bas2.params[1].evennum  ==  3
        @test_throws ErrorException CompEcon.Basis(CompEcon.Lin(), [0.0], 5)
        @test_throws ErrorException CompEcon.Basis(CompEcon.Lin(), [0.0,1.5,4.0], 5)
    end

    @testset "test constructor" begin

       valbas_lin = @inferred CompEcon.evalbase(bas1.params[1])
       valbas_spl = @inferred CompEcon.evalbase(CompEcon.SplineParams([0.0,1.0,2.0,3.0,4.0],0,1))

       @test  valbas_lin  ==  valbas_spl

    end

    @testset "test derivative" begin

    	#test derivative of basis indirectly, by fitting exponential function
		bas2 = Basis(LinParams([0,1.0],1000000))
		coeffs = funfitf(bas2,exp)
		points = rand(10000)

		d1 = funeval(coeffs,bas2,points,1)

        @test_approx_eq_eps d1 exp(points) 1e-7

    end

end
