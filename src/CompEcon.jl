VERSION >= v"0.4.0-dev+6521" && __precompile__()

module CompEcon

#=
TODO: Maybe package in two modules, CompEcon, and CompEcon.(Original|Classic)

TODO: still need to write fund, minterp
TODO: also need splidop, lindop
TODO: funeval fails for scalar input and does weird thing for 1-element
      vector input

TODO: remove the x return value for evalbase
TODO: Simplify API by only allowing `order` in evalbase to be Int. Then we get
type stability in the first (soon to be only) return argument. We could
implement another method that expects `order=Vector`

TODO: rethink the `order::Matrix` thing. Maybe a convention we can adopt is
the following:

- If calling evalbase on a `*Params` type, we will say that order must be of
type `Integer` or `T<:Integer, Vector{T}`. The implementation for the `Vector`
case will just loop over each element and then call the `Integer` method,
which will have all the code
- If calling `evalbase` on a `Basis{1}` type, the semantics will be the same
as for the `*Params` case
- If calling `evalbase` on a `Basis{N>1}` type, then `order` can either be a
`T<:Integer Vector{T}` or `T<:Integer Matrix{T}`. If it is a vector, we call
the `Int` evalbase for the `*Params` along each dimension, then combine as
necessary. If it is a `Matrix` then we say that each column is the orders that
should be computed for each dimension and we farm it out that way.

So method signatures would be

```
evalbase(p::ChebParam, order::Int) => Matrix{Float64}
evalbase(p::SplineParam, order::Int) => SparseMatrixCSC{Float64}
evalbase(p::LinParam, order::Int) => SparseMatrixCSC{Float64}

evalbase(p::ChebParam, order::Vector{Int}) => Vector{Matrix{Float64}}
evalbase(p::SplineParam, order::Vector{Int}) => Vector{SparseMatrixCSC{Float64}}
evalbase(p::LinParam, order::Vector{Int}) => Vector{SparseMatrixCSC{Float64}}

evalbase(b::Basis{1}, order::Vector{Int}) => .... can't remember
evalbase(b::Basis{N}, order::Matrix{Int}) => Vector{... can't remembers}
```

TODO: potentially add another abstract type
      `abstract AbstractSpline <: BasisFamily` and then make
      `Lin <: AbstractSpline` and `Spline <: AbstractSpline` So I can preserve
      the type info for the Base.SparseMatrix.SparseMatrixCSC{Float64,Int64}
      TM type parameter on the `BasisStructure` for a hybrid Spline/Lin Basis

=#


#=
Note that each subtype of `BT<:BasisFamily` (with associated `PT<:BasisParam`)
will define the following constructor methods:

```julia
# basis constructors
Basis(::BT, args...)
Basis(::PT)

# node constructor
nodes(::PT)
```

=#

using Compat

import Base: ==

export golden_method

# zeros
export bisect, brenth, brent, ridder, expand_bracket, divide_bracket

# types
export BasisFamily, Cheb, Lin, Spline, Basis,
       BasisParams, ChebParams, LinParams, SplineParams,
       AbstractBasisStructureRep, Tensor, Expanded, Direct,
       BasisStructure, Interpoland

# functions
export old_name, nodes, revert, get_coefs, funfitxy, funfitf, funeval,
       derivative_op, funbasex, row_kron, evaluate, fit!, update_coefs!,
       complete_polynomial, complete_polynomial!


# old API only
export fundef, fundefn, funnode, funbase, funbasex, funeval, funbconv,
       chebdef, chebnode, chebbase, chebbasex, chebdop,
       splidef, splinode, splibase, splibasex, splidop,
       lindef, linnode, linbase, lindop

# quad
export qnwlege, qnwcheb, qnwsimp, qnwtrap, qnwbeta, qnwgamma, qnwequi, qnwnorm,
       qnwunif, qnwlogn, qnwgh, qnwmonomial,
       quadrect,
       gridmake,
       do_quad

include("util.jl")
include("zeros.jl")
include("optimization.jl")
include("original/core.jl")
include("basis.jl")
include("basis_structure.jl")
include("interp.jl")

# include the rest of the original API
include("ml_compat.jl")
include("original/cheb.jl")
include("original/spli.jl")
include("original/lin.jl")

# include the rest of the Julian API
include("cheb.jl")
include("spline.jl")
include("lin.jl")
include("complete.jl")

# include quad
include("quad.jl")


end # module
