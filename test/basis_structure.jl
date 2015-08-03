# Tests conversion among basis-structure representations
module TestBasis_Structure

using CompEcon
using Base.Test
using FactCheck


facts("Test Basis Structure Representations") do
    bt = [Spline(), Cheb()]
    n = [9, 7]
    a = [-3, 1e-4]
    b = [3, 10.0]
    params = [SplineParams(linspace(-3, 3, 9), 0, 1), ChebParams(7, 1e-4, 10.0)]

    # construct multivariate basis
    mb = Basis{2,BasisFamily,BasisParams}(bt, n, a, b, params)

    # construct evaluation points
    x1 = linspace(a[1],b[1],15)
    x2 = linspace(a[2],b[2],10)
    X = [kron(ones(10),x1) kron(x2,ones(15))]

    # construct direct basis-structure representation
    PHIdirect = BasisStructure(mb,Direct(),X)

    # construct tensor basis-structure representation
    PHIdirect = BasisStructure(mb,Tensor(),X)

    # construct direct basis-structure representation
    PHIdirect = BasisStructure(mb,Expanded(),X)


    convert(bst::Type{Expanded}, bs::BasisStructure{Direct})

    context("constructors") do
        @fact map(x->isa(x, Basis{1}), (b1, b2)) --> (true, true)
        @fact isa(b_both, Basis{2}) --> true

        # use convenience outer constructor to not supply type parameter when
        # constructing multi-varaite basis
        @fact Basis(bt, n, a, b, params) --> b_both

        # use Basis(Type, args...) constructor
        @fact Basis(Spline, CompEcon.old_params(params[1])...) --> b1
        @fact Basis(Cheb, CompEcon.old_params(params[2])...) --> b2

        # use Basis(params) constructor
        @fact Basis(params[1]) --> b1
        @fact Basis(params[2]) --> b2

        # test Basis(b1, b2) format
        @fact Basis(b1, b2) --> b_both

        # test unicode ×
        @fact b1 × b2 --> b_both

        # test Basis(params1, params2)
        @fact Basis(params...) --> b_both

        # test fundefn type methods
        @fact b_spline2d --> Basis(b1, b1)
        @fact b_cheb2d --> Basis(b2, b2)

        # test that basis of different dimensions are not equal
        @fact ==(b_spline2d, b1) --> false

        # test concrete types come out when possible
        @fact isa(Basis(b1, b1), Basis{2,Spline,SplineParams}) --> true
        @fact isa(Basis(b2, b2), Basis{2,Cheb,ChebParams}) --> true
        @fact isa(Basis(b2, b2, b1), Basis{3,BasisFamily,BasisParams}) --> true
    end