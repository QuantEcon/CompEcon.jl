# Tests constructors and subtype relationships
module TestBasis

using CompEcon
using Base.Test
using FactCheck


facts("Test Basis") do
    bt = [Spline(), Cheb()]
    n = [9, 7]
    a = [-3, 1e-4]
    b = [3, 10.0]
    params = [SplineParams(linspace(-3, 3, 9), 0, 1), ChebParams(7, 1e-4, 10.0)]

    # use univariate constructors
    b1, b2 = map(Basis, bt, n, a, b, params)

    # directly construct multivariate basis using (default) inner constructor
    b_both = Basis{2}(bt, n, a, b, params)

    context("constructors") do
        @fact map(x->isa(x, Basis{1}), (b1, b2)) => (true, true)
        @fact isa(b_both, Basis{2}) => true

        # use convenience outer constructor to not supply type parameter when
        # constructing multi-varaite basis
        @fact Basis(bt, n, a, b, params) => b_both

        # use Basis(Type, args...) constructor
        @fact Basis(Spline, CompEcon.old_params(params[1])...) => b1
        @fact Basis(Cheb, CompEcon.old_params(params[2])...) => b2

        # use Basis(params) constructor
        @fact Basis(params[1]) => b1
        @fact Basis(params[2]) => b2

        # test Basis(b1, b2) format
        @fact Basis(b1, b2) => b_both

        # test unicode ×
        @fact b1 × b2 => b_both

        # test Basis(params1, params2)
        @fact Basis(params...) => b_both
    end

    context("getindex and combining preserves basis") do
        @fact b_both[1] => b1
        @fact b_both[2] => b2
        @fact Basis(b_both[1], b_both[2]) => b_both
        for b in [b1, b2]
            @fact b[1] => b
            @fact_throws ErrorException b[2]
        end
    end

    context("ndims, length, and size") do
        for b in [b1, b2]
            @fact ndims(b) => 1
            @fact length(b) => b.n[1]
            @fact size(b) => (b.n[1],)
            @fact size(b, 1) => b.n[1]
            @fact_throws ErrorException size(b, 2)
        end

        @fact ndims(b_both) => 2
        @fact length(b_both) => prod(n)
        @fact size(b_both) => tuple(n...)
        @fact size(b_both, 1) => n[1]

    end

    context("test nodes") do
        # extract nodes from individual basis
        n1 = nodes(b1)[1]
        n2 = nodes(b2)[1]

        # extract product nodes as well as nodes along both dimensions using
        # the 2d basis
        n_both, (n1_both, n2_both) = nodes(b_both)

        # test the nodes from basis 1 are what we expect, are the same as the
        # corresponding nodes from the 2d basis and have the correct length
        @fact n1 => collect(linspace(-3, 3, 9))
        @fact n1 => n1_both
        @fact length(n1) => b1.n[1]

        # test that the nodes from basis 2 are the same as corresponding nodes
        # from 2d basis and have correct length
        @fact n2 => n2_both
        @fact length(n2) => b2.n[1]

        # verify that the nodes from combined 2d basis is correct size
        @fact size(n_both) => (length(n1)*length(n2), 2)
    end
end

end
