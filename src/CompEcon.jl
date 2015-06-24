module CompEcon

#=
TODO: Maybe package in two modules, CompEcon, and CompEcon.(Original|Classic)

TODO: still need to write fund, minterp
TODO: also need splidop, lindop
TODO: funeval fails for scalar input and does weird thing for 1-element
      vector input
TODO: Decide if I should move core algorithms into julian api files or leave
      them in src/original. Right now src/original/core repeats a lot of code
      that is better written in src/
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

export golden_method

# types
export BasisFamily, Cheb, Lin, Spline, Basis,
       BasisParams, ChebParams, LinParams, SplineParams,
       AbstractBasisStructureRep, Tensor, Expanded, Direct,
       BasisStructure

# functions
export old_name, nodes, revert, get_coefs, funfitxy, funfitf, funeval,
       derivative_op, funbasex, row_kron


# old API only
export fundef, fundefn, funnode, funbase, funbasex, funeval, funbconv,
       chebdef, chebnode, chebbase, chebbasex, chebdop,
       splidef, splinode, splibase, splibasex, splidop,
       lindef, linnode, linbase, lindop

include("util.jl")
include("optimization.jl")
include("original/core.jl")
include("basis.jl")
include("basis_structure.jl")
include("interp.jl")

# include the rest of the original API
include("original/cheb.jl")
include("original/spli.jl")
include("original/lin.jl")

# include the rest of the Julian API
include("cheb.jl")
include("spline.jl")
include("lin.jl")


end # module
