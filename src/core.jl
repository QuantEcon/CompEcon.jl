# TODO: When all is said and done, change `parms` to `params`
# TODO: Maybe package in two modules, CompEcon, and CompEcon.(Original|Classic)
# TODO: Write a test that has two Basis objects, separates them with getindex
#       then puts them back together with a constructor or × and see if we get
#       the same thing out

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

old_name(::Lin) = :lin
old_name(::Cheb) = :cheb
old_name(::Spline) = :spli

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

old_parms(p::LinParms) = Any[p.breaks, p.evennum]
old_parms(p::ChebParms) = Any[p.n, p.a, p.b]
old_parms(p::SplineParms) = Any[p.breaks, p.evennum, p.k]

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

# constructor to allow `Basis(Spline, breaks, evennum, k)` instead of just
# `Basis(Spline(), breaks, evennum, k)`
Basis{BT<:BasisFamily}(::Type{BT}, args...) = Basis(BT(), args...)

# combining basis
function Basis(b1::Basis, bs::Basis...)  # fundef-esque method
    Nb = length(bs)
    basistype = vcat(b1.basistype, [bs[i].basistype for i=1:Nb]...)
    n = vcat(b1.n, [bs[i].n for i=1:Nb]...)
    b = vcat(b1.b, [bs[i].b for i=1:Nb]...)
    a = vcat(b1.a, [bs[i].a for i=1:Nb]...)
    parms = vcat(b1.parms, [bs[i].parms for i=1:Nb]...)
    N = length(n)
    Basis{N}(basistype, n, a, b, parms)
end

Base.×(b1::Basis, b2::Basis) = Basis(b1, b2)

# construct basis out of multiple Params
Basis(p::AnyParm, ps::AnyParm...) = Basis(Basis(p),
                                          Basis[Basis(p) for p in ps]...)


# convert old API to new API
function Base.convert(::Type{Basis}, b::Dict{Symbol, Any})
    Basis{b[:d]}(b[:basetype], b[:n], b[:a], b[:b], b[:parms])
end

# convert new API to old API
function revert(b::Basis)
    B = Dict{Symbol, Any}()
    B[:d] = ndims(b)
    B[:n] = b.n
    B[:a] = b.a
    B[:b] = b.b
    B[:basetype] = Symbol[old_name(bt) for bt in b.basistype]
    B[:parms] = Any[old_parms(p) for p in b.parms]
    B
end

# fundefn methods
Basis(bt::BasisFamily, n::Matrix, a::Matrix, b::Matrix, order::Int=3) =
    Basis(fundefn(old_name(bt), n, a, b, order))

Basis(bt::BasisFamily, n::Vector, a::Vector, b::Vector, order::Int=3) =
    Basis(bt, n[:, :], a[:, :], b[:, :], order)


# separating Basis
function Base.getindex{N}(basis::Basis{N}, n::Int)
    n < 0 || n > N && error("n must be between 1 and $N")
    return Basis{1}(basis.basistype[[n]],  # double `[[` to retain Vector
                    basis.n[[n]],
                    basis.b[[n]],
                    basis.a[[n]],
                    basis.parms[[n]])
end

# Define standard Julia methods for Basis
Base.ndims{N}(::Basis{N}) = N
Base.length(b::Basis) = prod(b.n)
Base.size(b::Basis, n::Int) = length(b[n])
Base.size{N}(b::Basis{N}) = tuple(b.n...)::NTuple{N, Int64}

# ------------------- #
# BasisStructure Type #
# ------------------- #

abstract AbstractBasisStructureRep
immutable Tensor <: AbstractBasisStructureRep end
immutable Direct <: AbstractBasisStructureRep end
immutable Expanded <: AbstractBasisStructureRep end

typealias ABSR AbstractBasisStructureRep

immutable BasisStructure{BST<:AbstractBasisStructureRep}
    order::Matrix{Int}
    vals::Array{AbstractMatrix}
end

Base.ndims(bs::BasisStructure) = size(bs.order, 2)
# TODO: determine if I want `size(bs, i::Int)` methods to do
# `prod([size(b.vals[1,j], i) for j=1:ndims(b)])`

# common checks for all convert methods
function check_convert{BST<:ABSR}(::Type{BST}, bs::BasisStructure, order)
    d = ndims(bs)
    numbas, d1 = size(bs.order)

    d1 != d && error("ORDER incompatible with basis functions")  # 35-37

    # 39-41
    if any(minimum(order) .< bs.order)
        error("Order of derivative operators exceeds basis")
    end
    return d, numbas, d1
end

Base.convert(bst::Type{Expanded}, bs::BasisStructure{Direct}, args...) = bs


# funbconv from direct to expanded
function Base.convert(bst::Type{Expanded}, bs::BasisStructure{Direct},
                      order=fill(0, 1, size(bs.order, 2)))
    d, numbas, d1 = check_convert(bst, bs, order)
    n = prod([size(bs.vals[1, j], 2) for j=1:d])

    vals = Array(AbstractMatrix, numbas)

    for i=1:numbas
        vals[i] = bs.vals[order[i, d] - bs.order[d]+1, d]  # 63
        for j=d-1:-1:1
            vals[i] = row_kron(vals[i],
                                   bs.vals[order[i, j] - bs.order[j]+1, j])  #65
        end
    end
    BasisStructure{Expanded}(order, vals)
end

# funbconv from tensor to expanded
function Base.convert(bst::Type{Expanded}, bs::BasisStructure{Tensor},
                      order=fill(0, 1, size(bs.order, 2)))
    d, numbas, d1 = check_convert(bst, bs, order)
    m = prod([size(bs.vals[1, j], 1) for j=1:d])
    n = prod([size(bs.vals[1, j], 2) for j=1:d])

    vals = Array(AbstractMatrix, numbas)

    for i=1:numbas  # 54
        vals[i] = bs.vals[order[i, d] - bs.order[d]+1, d]  # 55
        for j=d-1:-1:1  # 56
            # 57
            vals[i] = kron(vals[i], bs.vals[order[i, j] - bs.order[j]+1, j])
        end
    end

    BasisStructure{Expanded}(order, vals)
end

# funbconv from tensor to direct
function Base.convert(bst::Type{Direct}, bs::BasisStructure{Tensor},
                      order=fill(0, 1, size(bs.order), 2))
    d, numbas, d1 = check_convert(bst, bs, order)
    vals = Array(AbstractMatrix, numbas, d)
    ind = cell(1, d)  # 78

    for j=1:d, i=1:size(bs.vals, 1)
        if !(isempty(bs.vals[i, j]))  # 81
            ind[j] = collect(1:size(bs.vals[i, j], 1))  # 82
            break
        end
    end

    ind = gridmake(ind...)  # 87

    for j=1:d, i=1:numbas
        if !(isempty(bs.vals[i, j]))  # 90
            vals[i, j] = vals[i, j][ind[j], :]
        end
    end

    BasisStructure{Direct}(order, vals)
end

# ------------------------------------- #
# Main functionality at the `fun` level #
# ------------------------------------- #

function nodes(b::Basis)  # funnode method
    d = ndims(b)
    xcoord = Vector{Float64}[nodes(b.parms[j]) for j in 1:d]
    x = gridmake(xcoord...)
    return x, xcoord
end


# code to be run at the top of each `BasisStructure` constructor
# it enforces compatibility of arguments and computes common items
function check_basis_structure(basis::Basis, x, order, bformat)
    d = ndims(basis)
    if d > 1 && size(order, 2) == 1  # 62
        order = order * ones(Int, 1, d)  # 63
    end  # 64

    if d == 1 && isa(order, Int)
        order = fill(order, 1, 1)
    end

    # initialize basis structure (66-74)
    m = size(order, 1)
    if m > 1
        minorder = fill(minimum(order), 1, d)
        numbases = (maximum(order) - minorder) + 1
    else
        minorder = order + zeros(Int, 1, d)
        numbases = fill(1, 1, d)
    end

    return d, m, order, minorder, numbases
end

# quick function to take order+vals and return expanded form for 1d problems
function to_expanded(out_order::Matrix{Int}, vals::Array)
    vals = vals[collect(out_order + (1 - minimum(out_order)))]
    BasisStructure{Expanded}(out_order, vals)
end

# method to construct BasisStructure in direct or expanded form based on
# a matrix of `x` values
function BasisStructure(basis::Basis,
                        x::Array{Float64}=nodes(basis)[1], order=0,
                        bformat::Direct=Direct())  # funbasex

    d, m, order, minorder, numbases = check_basis_structure(basis, x, order,
                                                            bformat)
    # 76-77
    out_order = minorder
    out_format = bformat
    vals = Array(AbstractMatrix, maximum(numbases), d)

    # now do direct form, will convert to expanded later if needed
    for j=1:d
        # 126-130
        if (m > 1)
            orderj = unique(order[:, j])
        else
            orderj = order[1, j]
        end

        #131-135
        if length(orderj) == 1
            vals[1, j] = evalbase(basis.parms[j], x[:, j], orderj)[1]
        else
            vals[orderj-minorder[j]+1, j] =
                evalbase(basis.parms[j], x[:, j], orderj)[1]
        end
    end

    # if d == 1, switch to expanded format and return it directly
    # 140-145
    if size(vals, 2) == 1  # 1 dimension
        return to_expanded(order, vals)
    end

    # construct Direct Format
    BasisStructure{Direct}(out_order, vals)
end

function BasisStructure(basis::Basis,
                        x::Array{Float64}=nodes(basis)[1], order=0,
                        bformat::Expanded=Expanded())  # funbasex
    # create direct form, then convert to expanded
    convert(Expanded, BasisStructure(basis, x, order, Direct()))
end

# TODO: figure out a better way than just x::Array{Any} to do this,
#       Maybe x::Union(Array{Vector}, Union(Any)) because nodes(b)[2]
#       will Array{Vector}. Seems foolish to throw away that type information
#       just so dispatch will reach this function
function BasisStructure(basis::Basis, x::Array{Any}, order=0,
                        bformat::Tensor=Tensor())  # funbasex
    d, m, order, minorder, numbases = check_basis_structure(basis, x, order,
                                                        bformat)
    out_order = minorder
    out_format = bformat
    vals = Array(AbstractMatrix, maximum(numbases), d)

    # construct tensor base
    for j=1:d
        # 113-117
        if (m > 1)
            orderj = unique(order[:, j])
        else
            orderj = order[1, j]
        end

        #118-122
        if length(orderj) == 1
            vals[1, j] = evalbase(basis.parms[j], x[j], orderj)[1]
        else
            vals[orderj-minorder[j]+1, j] = evalbase(basis[:parms][j], x[j],
                                                     orderj)[1]
        end
    end

    size(vals, 2) == 1 ? to_expanded(order, vals) :
                         BasisStructure{Tensor}(out_order, vals)
end

# add method to funbasex that calls the above constructors
funbasex(basis::Basis, x=nodes(basis)[1], order=0, bformat::ABSR=Direct()) =
    BasisStructure(basis, x, order, bformat)

funbase(basis::Basis, x=nodes(basis)[1], order=fill(0, 1, ndims(basis))) =
    funbasex(basis, x, order, Expanded()).vals[1]


function get_coefs(basis::Basis, bs::BasisStructure{Tensor}, y)
    if any(bs.order[1, :] != 0)
        error("invalid basis structure - first elements must be order 0")
    end
    bs.vals = bs.vals[1, :]  # 68
    ckronxi(bs.vals, y, ndims(basis):-1:1)  # 66
end


get_coefs(basis::Basis, bs::BasisStructure{Direct}, y) =
    get_coefs(basis, convert(Expanded, bs), y)

get_coefs(basis::Basis, bs::BasisStructure{Expanded}, y) = bs.vals[1] \ y

# common checks to be run at the top of each funfit
function check_funfit(basis::Basis, x, y)
    m = size(y, 1)
    length(basis) > m && error("Can't be more basis funcs than points in y")
    return m
end

function funfitxy(basis::Basis, bs::BasisStructure, y)
    m = check_funfit(basis, bs, y)
    c = get_coefs(basis, bs, y)
    c, bs
end

function funfitxy(basis::Basis, x::Array{Any}, y)
    m = check_funfit(basis, x, y)

    # additional checks for cell array
    mm = prod([size(x[i], 1) for i=1:size(x, 2)])
    mm != m && error("x and y are incompatible")

    # get Tensor form
    bs = BasisStructure(basis, x, 0)

    # get coefs and return
    c = get_coefs(basis, bs, y)
    c, bs
end

function funfitxy(basis::Basis, x, y)
    # check input sizes
    m = check_funfit(basis, x, y)

    # additional check
    size(x, 1) != m && error("x and y are incompatible")

    # compute expanded basis structure, get cofs, and return
    bs = BasisStructure(basis, x, 0, Expanded())
    c = get_coefs(basis, bs, y)
    c, bs
end

function funfitf(basis::Basis, f::Function, args...)
    x = nodes(basis)[1]
    y = f(x, args...)
    funfitxy(basis, x, y)[1]
end

# funeval wants to evaluate at a matrix. As a stop-gap until I find some
# time, this method makes a scalar x into a 1x1 matrix
funeval(c, basis::Basis, x::Real, order=0) =
    funeval(c, basis, fill(x, 1, 1), order)

# similar to above for vectors (size will be nx1)
funeval(c, basis::Basis, x::Vector, order=0) =
    funeval(c, basis, x[:, :], order)

function funeval(c, basis::Basis, x::Matrix, order=0)
    d = ndims(basis)
    if size(x, 2) != d
        error("x must have d=$(d) columns")  # 62
    end

    if size(order, 2) == 1
        order = order * fill(1, 1, d)  # 65
    end
    B = BasisStructure(basis, x, order)  # 67

    funeval(c, B, order)
end

function funeval(c, bs::BasisStructure{Tensor},
                 order::Matrix{Int}=fill(0, 1, size(bs.order, 2)))  # funeval1
    kk, d = size(order)  # 95
    # 98 reverse the order of evaluation: bs(d)xB(d-1)x...xB(1)
    order = flipdim(order+1*(size(bs.vals, 1) * (0:d-1)' - bs.order+1), 2)

    # 99
    nx = prod([size(bs.vals[1, j], 1) for j=1:d])

    f = zeros(nx, size(c, 2), kk)  # 100

    for i=1:kk
        f[:, :, i] = ckronx(bs.vals, c, order[i, :])  # 102
    end
    return squeeze_trail(f)
end

function funeval(c, bs::BasisStructure{Direct},
                 order::Matrix{Int}=fill(0, 1, size(bs.order, 2)))  # funeval2
    kk, d = size(order)  # 95
    # 114 reverse the order of evaluation: B(d)xB(d-1)x...xB(1)
    order = flipdim(order+1*(size(bs.vals, 1) * (0:d-1)' - bs.order+1), 2)

    f = zeros(size(bs.vals[1], 1), size(c, 2), kk)  # 116

    for i=1:kk
        f[:, :, i] = cdprodx(bs.vals, c, order[i, :])  # 118
    end
    return squeeze_trail(f)
end

# TODO: bugs here. Need to find a better way to replicate the `iscell`
#       calls in Matlab because bs.vals is never and Array{Any}
function funeval(c, bs::BasisStructure{Expanded},
                 order::Matrix{Int}=fill(0, 1, size(bs.order, 2)))  # funeval3

    if isempty(order)
        # TODO: this is dead code because b.vals <: Array{AbstractMatrix}
        if isa(bs.vals, Array{Any})
            kk = length(bs.vals)  # 126
            order = 1:kk'  # 127
        else
            kk = 1  # 128
            order = 1  # 129
        end
    else
        kk = size(order, 1)  # 133
    end

    if isa(bs.vals, Array{Any})
        nx = size(bs.vals[1], 1)
        f = zeros(nx, size(c, 2), kk)
        for i=1:kk
            # 140 determine which element of bs.vals is the desired basis
            ii = Int[]
            for row=1:size(bs.order, 1)
                r = bs.order[row, :]
                if r == order[i, :]
                    push!(ii, row)
                end
            end

            # 141-143
            isempty(ii)  && error("Requested basis matrix is not available")

            length(ii) > 1 &&  warn("redundant request in funeval3")  # 145

            # NOTE: must do even when length[i] == 1 b/c want element of cell
            #       and indexing cell with vector in julia gives cell instead
            #       of the element
            ii = ii[1]  # 146

            f[:, :, i] = bs.vals[ii]*c  #148
         end
     else
        nx = size(bs.vals, 1)  # 151
        f = zeros(nx, size(c, 2), kk)  # 152
        for i=1:kk
            f[:, :, i] = bs.vals[1]*c  # 154
        end
    end

    return squeeze_trail(f)
end
