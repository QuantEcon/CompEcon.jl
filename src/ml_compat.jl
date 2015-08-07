# Stuff We need to maintain compatibility with the matlab api

## BasisFamily stuff

# needed for the `convert(::Basis, ::Dict)` method below
Base.convert(::Type{BasisFamily}, s::Symbol) =
    s == :spli ? Spline() :
    s == :cheb ? Cheb() :
    s == :lin  ? Lin() :
    error("Unknown basis type")

old_name(::Lin) = :lin
old_name(::Cheb) = :cheb
old_name(::Spline) = :spli

## BasisParams stuff

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


## Basis stuff

# convert old API to new API
function Base.convert(::Type{Basis}, b::Dict{Symbol, Any})
    btype = map(x->convert(BasisFamily, x), b[:basetype])
    param = BasisParams[convert(BasisParams, x) for x in b[:params]]
    Basis(btype, b[:n], b[:a], b[:b], param)
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


# add method to funbasex that creates a BasisStructure
funbasex(basis::Basis, x=nodes(basis)[1], order=0, bformat::ABSR=Direct()) =
    BasisStructure(basis, x, order, bformat)

funbase(basis::Basis, x=nodes(basis)[1], order=fill(0, 1, ndims(basis))) =
    funbasex(basis, x, order, Expanded()).vals[1]
