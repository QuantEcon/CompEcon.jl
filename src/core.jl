module AA

# TODO: When all is said and done, change `parms` to `params`
# TODO: Maybe package in two modules, CompEcon, and CompEcon.(Original|Classic)
# TODO: Write a test that has two Basis objects, separates them with getindex
#       then puts them back together with a constructor or × and see if we get
#       the same thing out

include("util.jl")
include("original/core.jl")

# -------- #
# Families #
# -------- #

abstract BasisFamily
immutable Cheb <: BasisFamily end
immutable Lin <: BasisFamily end
immutable Spline <: BasisFamily end


typealias AnyFamily Union(Cheb, Lin, Spline)

Base.convert(::Type{AnyFamily}, s::Symbol) =
    s == :spli ? Spline() :
    s == :cheb ? Cheb() :
    s == :lin  ? Lin() :
    error("Unknown basis type")

# ---------- #
# Parameters #
# ---------- #

abstract BasisParms
immutable ChebParms <: BasisParms
    n::Int
    a::Int
    b::Int
end

immutable SplineParms <: BasisParms
    breaks::Vector{Float64}
    evennum::Int
    k::Int
end

immutable LinParms <: BasisParms
    breaks::Vector{Float64}
    evennum::Int
end

typealias AnyParm Union(ChebParms, SplineParms, LinParms)

Base.convert(::Type{AnyParm}, p::Vector{Any}) =
    length(p) == 3 && isa(p[1], Number) ? ChebParms(p...) :
    length(p) == 3 ? SplineParms(p...) :
    length(p) == 2 ? LinParms(p...) :
    error("Unknown parameter type")

old_name(::LinParms) = :lin
old_name(::ChebParms) = :cheb
old_name(::SplineParms) = :spli
old_name(::Lin) = :lin
old_name(::Cheb) = :cheb
old_name(::Spline) = :spli

# ---------- #
# Basis Type #
# ---------- #

immutable Basis{N}
    basistype::Vector{AnyFamily}
    n::Vector{Int}
    a::Vector{Float64}
    b::Vector{Float64}
    parms::Vector{AnyParm}
end

# univariate basis. Helper to wrap all args in arrays
Basis(bt::BasisFamily, n::Int, a::Float64, b::Float64, parms::BasisParms) =
    Basis{1}([bt], [n], [a], [b], [parms])

# combining basis
function Basis(bs::Basis...)
    Nb = length(bs)
    basistype = vcat([bs[i].basistype for i=1:Nb]...)
    n = vcat([bs[i].n for i=1:Nb]...)
    b = vcat([bs[i].b for i=1:Nb]...)
    a = vcat([bs[i].a for i=1:Nb]...)
    parms = vcat([bs[i].parms for i=1:Nb]...)
    N = length(n)
    Basis{N}(basistype, n, a, b, parms)
end

Base.×(b1::Basis, b2::Basis) = Basis(b1, b2)

function Base.convert(::Type{Basis}, b::Dict{Symbol, Any})
    Basis{b[:d]}(b[:basetype], b[:n], b[:a], b[:b], b[:parms])
end

# fundefn
Basis(bt::BasisFamily, n::Matrix, a::Matrix, b::Matrix, order::Int=3) =
    Basis(fundefn(old_name(bt), n, a, b, order))

Basis(bt::BasisFamily, n::Vector, a::Vector, b::Vector, order::Int=3) =
    Basis(bt, n[:, :], a[:, :], b[:, :], order)


Base.ndims{N}(basis::Basis{N}) = N

# separating Basis
function Base.getindex{N}(basis::Basis{N}, n::Int)
    n < 0 || n > N && error("n must be between 1 and $N")
    return Basis{1}(basis.basistype[[n]],  # double `[[` to retain Vector
                    basis.n[[n]],
                    basis.b[[n]],
                    basis.a[[n]],
                    basis.parms[[n]])
end

function nodes(b::Basis)
    d = ndims(b)
    xcoord = Vector[nodes(b.parms[j]) for j in 1:d]
    x = gridmake(xcoord...)
    return x, xcoord
end

include("spline.jl")
include("cheb.jl")

end  # module
