# -------- #
# Families #
# -------- #

abstract BasisFamily
immutable Cheb <: BasisFamily end
immutable Lin <: BasisFamily end
immutable Spline <: BasisFamily end

typealias AnyFamily Union(Cheb, Lin, Spline)

# needed for the `convert(::Basis, ::Dict)` method below
Base.convert(::Type{BasisFamily}, s::Symbol) =
    s == :spli ? Spline() :
    s == :cheb ? Cheb() :
    s == :lin  ? Lin() :
    error("Unknown basis type")

old_name(::Lin) = :lin
old_name(::Cheb) = :cheb
old_name(::Spline) = :spli

# ---------- #
# Parameters #
# ---------- #

abstract BasisParams
immutable ChebParams <: BasisParams
    n::Int
    a::Float64
    b::Float64
end

function Base.writemime(io::IO, ::MIME"text/plain", p::ChebParams)
    m = string("Chebyshev interpoland parameters with ",
               "$(p.n) basis functions from $(p.a), $(p.b)")
    print(m)
end

immutable SplineParams <: BasisParams
    breaks::Vector{Float64}
    evennum::Int
    k::Int
end

function Base.writemime(io::IO, ::MIME"text/plain", p::SplineParams)
    m = string("$(p.k) order spline interpoland parameters from ",
               "$(p.breaks[1]), $(p.breaks[end])")
    print(m)
end

immutable LinParams <: BasisParams
    breaks::Vector{Float64}
    evennum::Int
end

function Base.writemime(io::IO, ::MIME"text/plain", p::LinParams)
    m = string("$(p.k) order piecewise linear interpoland parameters ",
               "from $(p.breaks[1]), $(p.breaks[end])")
    print(m)
end

typealias AnyParam Union(ChebParams, SplineParams, LinParams)

# needed for the `convert(::Basis, ::Dict)` method below
Base.convert(::Type{BasisParams}, p::Vector{Any}) =
    length(p) == 3 && isa(p[1], Number) ? ChebParams(p...) :
    length(p) == 3 ? SplineParams(p...) :
    length(p) == 2 ? LinParams(p...) :
    error("Unknown parameter type")

old_name(::LinParams) = :lin
old_name(::ChebParams) = :cheb
old_name(::SplineParams) = :spli

old_params(p::LinParams) = Any[p.breaks, p.evennum]
old_params(p::ChebParams) = Any[p.n, p.a, p.b]
old_params(p::SplineParams) = Any[p.breaks, p.evennum, p.k]

# ---------- #
# Basis Type #
# ---------- #

immutable Basis{N}
    basistype::Vector{BasisFamily}  # Basis family
    n::Vector{Int}                  # number of points and/or basis functions
    a::Vector{Float64}              # lower bound of domain
    b::Vector{Float64}              # upper bound of domain
    params::Vector{BasisParams}     # params to construct basis
end

function Base.writemime{N}(io::IO, ::MIME"text/plain", b::Basis{N})
    m = """
    $N dimensional Basis on the hypercube formed by $(b.a) × $(b.b).
    Basis families are $(join(map(x->string(typeof(x)), b.basistype), " × "))
    """
    print(m)
end

# constructor that takes all arguments and ensures each has N elemnets
function Basis(basistype::Vector{BasisFamily}, n::Vector{Int},
               a::Vector{Float64}, b::Vector{Float64},
               params::Vector{BasisParams})
    N = length(basistype)
    @assert all(map(length, Any[n, a, b, params]) .== N)
    Basis{N}(basistype, n, a, b, params)
end

# univariate basis. Helper to wrap all args in arrays
Basis(bt::BasisFamily, n::Int, a::Float64, b::Float64, params::BasisParams) =
    Basis{1}([bt], [n], [a], [b], [params])

# constructor to allow `Basis(Spline, breaks, evennum, k)` instead of just
# `Basis(Spline(), breaks, evennum, k)`
Basis{BT<:BasisFamily}(::Type{BT}, args...) = Basis(BT(), args...)

# combining basis
function Basis(b1::Basis, bs::Basis...)  # fundef-esque method
    Nb = length(bs)
    basistype = vcat(b1.basistype, [bs[i].basistype for i=1:Nb]...)
    n = vcat(b1.n, [bs[i].n for i=1:Nb]...)
    a = vcat(b1.a, [bs[i].a for i=1:Nb]...)
    b = vcat(b1.b, [bs[i].b for i=1:Nb]...)
    params = vcat(b1.params, [bs[i].params for i=1:Nb]...)
    N = length(n)
    Basis{N}(basistype, n, a, b, params)
end

Base.×(b1::Basis, b2::Basis) = Basis(b1, b2)

# construct basis out of multiple Params
Basis(p::AnyParam, ps::AnyParam...) = Basis(Basis(p),
                                            Basis[Basis(p) for p in ps]...)


# convert old API to new API
function Base.convert(::Type{Basis}, b::Dict{Symbol, Any})
    Basis{b[:d]}(b[:basetype], b[:n], b[:a], b[:b], b[:params])
end

# convert new API to old API
function revert(b::Basis)
    B = Dict{Symbol, Any}()
    B[:d] = ndims(b)
    B[:n] = b.n
    B[:a] = b.a
    B[:b] = b.b
    B[:basetype] = Symbol[old_name(bt) for bt in b.basistype]
    B[:params] = Any[old_params(p) for p in b.params]
    B
end

# fundefn methods
Basis(bt::BasisFamily, n::Matrix, a::Matrix, b::Matrix, order::Int=3) =
    Basis(fundefn(old_name(bt), n, a, b, order))

# the [:, :] makes into Matrix and we call the previous method
Basis(bt::BasisFamily, n::Vector, a::Vector, b::Vector, order::Int=3) =
    Basis(bt, n[:, :], a[:, :], b[:, :], order)


# separating Basis
function Base.getindex{N}(basis::Basis{N}, n::Int)
    n < 0 || n > N && error("n must be between 1 and $N")
    return Basis{1}(basis.basistype[[n]],  # double `[[` to retain Vector
                    basis.n[[n]],
                    basis.a[[n]],
                    basis.b[[n]],
                    basis.params[[n]])
end

# Define standard Julia methods for Basis
Base.ndims{N}(::Basis{N}) = N
Base.length(b::Basis) = prod(b.n)
Base.size(b::Basis, i::Int) = length(b[i])  # uses method on previous line
Base.size{N}(b::Basis{N}) = tuple(b.n...)::NTuple{N, Int64}

# Bases of different dimension can't be equal
=={N,M}(::Basis{N}, ::Basis{M}) = false

# basis are equal if all fields of the basis are equal
=={N}(b1::Basis{N}, b2::Basis{N}) =
    all(map(nm->getfield(b1, nm) == getfield(b2, nm), fieldnames(b1)))

function nodes(b::Basis)  # funnode method
    d = ndims(b)
    xcoord = Vector{Float64}[nodes(b.params[j]) for j in 1:d]
    x = gridmake(xcoord...)
    return x, xcoord
end
