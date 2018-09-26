module CompEcon

using Reexport

using QuantEcon: gridmake, gridmake!, ckron, fix, fix!, qnwlege, qnwcheb,
                 qnwsimp, qnwtrap, qnwbeta, qnwgamma, qnwequi, qnwnorm,
                 qnwunif, qnwlogn, quadrect, do_quad

@reexport using BasisMatrices

import BasisMatrices: funeval, funfitxy, funfitf

# old API only
export fundef, fundefn, funnode, funbase, funbasex, funeval, funbconv,
       chebdef, chebnode, chebbase, chebbasex, chebdop,
       splidef, splinode, splibase, splibasex, splidop,
       lindef, linnode, linbase, lindop

# quad names
export qnwlege, qnwcheb, qnwsimp, qnwtrap, qnwbeta, qnwgamma, qnwequi, qnwnorm,
       qnwunif, qnwlogn, quadrect, do_quad

include("core.jl")
include("cheb.jl")
include("spli.jl")
include("lin.jl")
include("compat.jl")

end # module
