# Tests constructors and subtype relationships
module TestBasis

using CompEcon
using Base.Test
using FactCheck


facts("Test Basis") do
    bt = [Spline(), Cheb(), Lin()]
    n = [9, 7, 11]
    a = [-3, 1e-4, -4]
    b = [3, 10.0, 2]
    params = [SplineParams(linspace(-3, 3, 9), 0, 1), ChebParams(7, 1e-4, 10.0), LinParams(linspace(-4, 2, 11), 0)]

    # use univariate constructors
    b1, b2, b3 = map(Basis, bt, n, a, b, params)

    # directly construct multivariate basis using (default) inner constructor
    b_all = Basis{3,BasisFamily,BasisParams}(bt, n, a, b, params)

    # should be equivalent to Basis(b1, b1)
    b_spline3d = Basis(Spline(), [n[1], n[1]], [a[1], a[1]], [b[1], b[1]], [1, 1])

    # should be equivalent to Basis(b2, b2)
    b_cheb3d = Basis(Cheb(), [n[2], n[2]], [a[2], a[2]], [b[2], b[2]])

    # should be equivalent to Basis(b3, b3)
    b_lin3d = Basis(Lin(), [n[3], n[3]], [a[3], a[3]], [b[3], b[3]])

    context("constructors") do
        @fact map(x->isa(x, Basis{1}), (b1, b2, b3)) --> (true, true, true)
        @fact isa(b_all, Basis{3}) --> true

        # use convenience outer constructor to not supply type parameter when
        # constructing multi-varaite basis
        @fact Basis(bt, n, a, b, params) --> b_all

        # use Basis(Type, args...) constructor
        @fact Basis(Spline, Original.old_params(params[1])...) --> b1
        @fact Basis(Cheb, Original.old_params(params[2])...) --> b2
        @fact Basis(Lin, Original.old_params(params[3])...) --> b3

        # use Basis(params) constructor
        @fact Basis(params[1]) --> b1
        @fact Basis(params[2]) --> b2
        @fact Basis(params[3]) --> b3

        # test Basis(b1, b2, b3) format
        @fact Basis(b1, b2, b3) --> b_all

        # test unicode ×
        @fact b1 × b2 × b3 --> b_all

        # test Basis(params1, params2, params3)
        @fact Basis(params...) --> b_all

        # test fundefn type methods
        @fact b_spline3d --> Basis(b1, b1)
        @fact b_cheb3d --> Basis(b2, b2)
        @fact b_lin3d --> Basis(b3, b3)

        # test that basis of different dimensions are not equal
        @fact ==(b_spline3d, b1) --> false
        @fact ==(b_cheb3d, b2) --> false
        @fact ==(b_lin3d, b3) --> false

        # test concrete types come out when possible
        @fact isa(Basis(b1, b1), Basis{2,Spline,SplineParams}) --> true
        @fact isa(Basis(b2, b2), Basis{2,Cheb,ChebParams}) --> true
        @fact isa(Basis(b3, b3), Basis{2,Lin,LinParams}) --> true
        @fact isa(Basis(b2, b2, b1), Basis{3,BasisFamily,BasisParams}) --> true
        @fact isa(Basis(b3, b2, b1), Basis{3,BasisFamily,BasisParams}) --> true
    end

    context("getindex and combining preserves basis") do
        @fact b_all[1] --> b1
        @fact b_all[2] --> b2
        @fact b_all[3] --> b3
        @fact Basis(b_all[1], b_all[2], b_all[3]) --> b_all
        for b in [b1, b2, b3]
            @fact b[1] --> b
            @fact_throws ErrorException b[2]
        end
    end

    context("ndims, length, and size") do
        for b in [b1, b2, b3]
            @fact ndims(b) --> 1
            @fact length(b) --> b.n[1]
            @fact size(b) --> (b.n[1],)
            @fact size(b, 1) --> b.n[1]
            @fact_throws ErrorException size(b, 2)
        end

        @fact ndims(b_all) --> 3
        @fact length(b_all) --> prod(n)
        @fact size(b_all) --> tuple(n...)
        @fact size(b_all, 1) --> n[1]

    end

    context("test nodes") do
        # extract nodes from individual basis
        n1 = nodes(b1)[1]
        n2 = nodes(b2)[1]
        n3 = nodes(b3)[1]

        # extract product nodes as well as nodes along both dimensions using
        # the 3d basis
        n_all, (n1_all, n2_all, n3_all) = nodes(b_all)

        # test the nodes from basis 1 are what we expect, are the same as the
        # corresponding nodes from the 3d basis and have the correct length
        @fact n1 --> collect(linspace(-3, 3, 9))
        @fact n1 --> n1_all
        @fact length(n1) --> b1.n[1]

        # test that the nodes from basis 2 are the same as corresponding nodes
        # from 3d basis and have correct length
        @fact n2 --> n2_all
        @fact length(n2) --> b2.n[1]

        # test that the nodes from basis 3 are the same as corresponding nodes
        # from 3d basis and have correct length
        @fact n3 --> collect(linspace(-4, 2, 11))
        @fact length(n3) --> b3.n[1]


        # verify that the nodes from combined 3d basis is correct size
        @fact size(n_all) --> (length(n1)*length(n2)*length(n3), 3)
    end

    # call show (which calls writemime) so I can get 100% coverage :)
    context("Printing") do
        iob = IOBuffer())
        show(iob, ChebParams(10, -1, 1))
        show(iob, SplineParams(10, -1, 1))
        show(iob, LinParams(10, -1, 1))
        show(iob, b_all)
    end
end

end
