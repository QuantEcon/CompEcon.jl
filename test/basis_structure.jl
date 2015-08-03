# Tests conversion among basis-structure representations
module TestBasis_Structure

using CompEcon
using Base.Test
using FactCheck


facts("Test Basis Structure Representations") do
    mb = Basis(SplineParams(linspace(-3, 3, 9), 0, 1),
               ChebParams(7, 1e-4, 10.0))

    # construct evaluation points
    x1 = linspace(a[1],b[1],15)
    x2 = linspace(a[2],b[2],10)
    X = gridmake(x1, x2)

    # construct tensor basis-structure representation
    PHI_tensor = BasisStructure(mb, Tensor(), X)

    # construct direct basis-structure representation
    PHI_direct = BasisStructure(mb, Direct(), X)

    # construct direct basis-structure representation
    PHI_expanded = BasisStructure(mb, Expanded(), X)


    convert(bst::Type{Expanded}, bs::BasisStructure{Direct})

    context("constructors") do
    end

end
