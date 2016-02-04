# -------- #
# Families #
# -------- #

abstract BasisFamily
immutable Cheb <: BasisFamily end
immutable Lin <: BasisFamily end
immutable Spline <: BasisFamily end

# ---------- #
# Parameters #
# ---------- #

abstract BasisParams

type ChebParams <: BasisParams
    n::Int
    a::Float64
    b::Float64
end

function Base.show(io::IO, p::ChebParams)
    m = string("Chebyshev interpoland parameters with ",
               "$(p.n) basis functions from $(p.a), $(p.b)")
    print(io, m)
end

type SplineParams <: BasisParams
    breaks::Vector{Float64}
    evennum::Int
    k::Int
end

# constructor to take a, b, n and form linspace for breaks
SplineParams(n::Int, a::Real, b::Real, k::Int=3) =
    SplineParams(collect(linspace(a, b, n)), 0, k)

function Base.show(io::IO, p::SplineParams)
    m = string("$(p.k) order spline interpoland parameters from ",
               "$(p.breaks[1]), $(p.breaks[end])")
    print(io, m)
end

type LinParams <: BasisParams
    breaks::Vector{Float64}
    evennum::Int
end

# constructor to take a, b, n and form linspace for breaks
LinParams(n::Int, a::Real, b::Real) =
    LinParams(collect(linspace(a, b, n)), 0)

function Base.show(io::IO, p::LinParams)
    m = string("piecewise linear interpoland parameters ",
               "from $(p.breaks[1]), $(p.breaks[end])")
    print(io, m)
end

for (T, TP) in [(Cheb, ChebParams), (Lin, LinParams), (Spline, SplineParams)]
    @eval _param(::$(T)) = $TP
    @eval _param(::Type{$(T)}) = $TP
end

# params of same type are equal if all their fields are equal
=={T<:BasisParams}(p1::T, p2::T) =
    all(map(nm->getfield(p1, nm) == getfield(p2, nm), fieldnames(p1)))::Bool

# Bases of different dimension can't be equal
=={T1<:BasisParams,T2<:BasisParams}(::T1, ::T2) = false

# ---------- #
# Basis Type #
# ---------- #

type Basis{N, BF<:BasisFamily, BP<:BasisParams}
    basistype::Vector{BF}  # Basis family
    n::Vector{Int}         # number of points and/or basis functions
    a::Vector{Float64}     # lower bound of domain
    b::Vector{Float64}     # upper bound of domain
    params::Vector{BP}     # params to construct basis
end

function Base.show{N}(io::IO, b::Basis{N})
    m = """
    $N dimensional Basis on the hypercube formed by $(b.a) × $(b.b).
    Basis families are $(join(map(x->string(typeof(x)), b.basistype), " × "))
    """
    print(io, m)
end

# constructor that takes all arguments and ensures each has N elemnets
function Basis{BF<:BasisFamily, BP<:BasisParams}(basistype::Vector{BF},
                                                 n::Vector{Int},
                                                 a::Vector{Float64},
                                                 b::Vector{Float64},
                                                 params::Vector{BP})
    N = length(basistype)
    @assert all(map(length, Any[n, a, b, params]) .== N)
    Basis{N,BF,BP}(basistype, n, a, b, params)
end

# univariate basis. Helper to wrap all args in arrays
function Basis{BF<:BasisFamily, BP<:BasisParams}(bt::BF, n::Int, a::Float64,
                                                 b::Float64, params::BP)
    Basis{1,BF,BP}([bt], [n], [a], [b], [params])
end

# constructor to allow `Basis(Spline, breaks, evennum, k)` instead of just
# `Basis(Spline(), breaks, evennum, k)`
Basis{BT<:BasisFamily}(::Type{BT}, args...) = Basis(BT(), args...)

# combining basis -- fundef-esque method
function Basis(b1::Basis, bs::Basis...)
    Nb = length(bs)
    basistype = vcat(b1.basistype, [bs[i].basistype for i=1:Nb]...)
    n = vcat(b1.n, [bs[i].n for i=1:Nb]...)
    a = vcat(b1.a, [bs[i].a for i=1:Nb]...)
    b = vcat(b1.b, [bs[i].b for i=1:Nb]...)
    params = vcat(b1.params, [bs[i].params for i=1:Nb]...)
    N = length(n)

    # determine if BF and BP type are concrete or use fall-back abstracts
    BF1, BP1 = map(x->typeof(x[1]), Any[b1.basistype, b1.params])
    BF = all(map(x->isa(x, BF1), basistype[2:end])) ? BF1 : BasisFamily
    BP = all(map(x->isa(x, BP1), params[2:end])) ? BP1 : BasisParams
    Basis{N,BF,BP}(basistype, n, a, b, params)::Basis{N,BF,BP}
end

Base.×(b1::Basis, b2::Basis) = Basis(b1, b2)

# construct basis out of multiple Params (type assertion for stability)
Basis(p::BasisParams, ps::BasisParams...) =
    Basis(Basis(p), Basis[Basis(p) for p in ps]...)::Basis{length(ps) + 1}

# fundefn type method
Basis{T<:BasisFamily}(bt::T, n::Vector, a::Vector, b::Vector) =
    Basis(map(_param(T), n, a, b)...)::Basis{length(n)}

# special method for Spline that adds `k` argument
Basis(::Spline, n::Vector, a::Vector, b::Vector, k::Vector=ones(Int, length(n))) =
    Basis(map(SplineParams, n, a, b, k)...)::Basis{length(n)}

# separating Basis -- just re construct it from the nth set of params
function Base.getindex{N}(basis::Basis{N}, n::Int)
    n < 0 || n > N && error("n must be between 1 and $N")
    Basis(basis.params[n])::Basis{1}
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
    all(map(nm->getfield(b1, nm) == getfield(b2, nm), fieldnames(b1)))::Bool

function nodes(b::Basis{1})
    x = nodes(b.params[1])
    (x, Vector{Float64}[x])
end

function nodes{N}(b::Basis{N})  # funnode method
    xcoord = Vector{Float64}[nodes(b.params[j]) for j in 1:N]
    x = gridmake(xcoord...)::Matrix{Float64}
    return x, xcoord
end
