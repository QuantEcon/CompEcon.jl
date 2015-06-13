module CompEcon

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
include("core.jl")

# include the rest of the original API
include("original/cheb.jl")
include("original/spli.jl")
include("original/lin.jl")

# include the rest of the Julian API
include("cheb.jl")
include("spline.jl")
include("lin.jl")


end # module
