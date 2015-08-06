# ------------------- #
# BasisStructure Type #
# ------------------- #

abstract AbstractBasisStructureRep
typealias ABSR AbstractBasisStructureRep

immutable Tensor <: ABSR end
immutable Direct <: ABSR end
immutable Expanded <: ABSR end

type BasisStructure{BST<:ABSR, TM<:AbstractMatrix}
    order::Matrix{Int}
    vals::Array{TM}
end

Base.writemime{BST}(io::IO, ::MIME"text/plain", b::BasisStructure{BST}) =
    print("BasisStructure{$BST} of order $(b.order)")

Base.ndims(bs::BasisStructure) = size(bs.order, 2)
# TODO: determine if I want `size(bs, i::Int)` methods to do
# `prod([size(b.vals[1,j], i) for j=1:ndims(b)])`

# not the same if either type parameter is different
function =={BST1<:ABSR,BST2<:ABSR}(::BasisStructure{BST1}, ::BasisStructure{BST2})
    false
end

function =={BST<:ABSR,TM1<:AbstractMatrix,TM2<:AbstractMatrix}(::BasisStructure{BST,TM1},
                                                               ::BasisStructure{BST,TM2})
    false
end

# if type parameters are the same, then it is the same if all fields are the
# same
function =={BST<:ABSR,TM<:AbstractMatrix}(b1::BasisStructure{BST,TM},
                                          b2::BasisStructure{BST,TM})
    b1.order == b2.order && b1.vals == b2.vals
end

# -------------- #
# Internal Tools #
# -------------- #

@inline function _checkx(N, x::AbstractMatrix)
    size(x, 2) != N && error("Basis is $N dimensional, x must have $N columns")
    x
end

@inline function _checkx{T<:Number}(N, x::AbstractVector{T})
    # if we have a 1d basis, we can evaluate at each point
    if N == 1
        return x
    end

    # If Basis is > 1d, one evaluation point and reshape to (1,N) if possible...
    if length(x) == N
        return reshape(x, 1, N)
    end

    # ... or throw an error
    error("Basis is $N dimensional, x must have $N elements")
end

@inline function _checkx{T<:Number}(N, x::Vector{Vector{T}})
    # for BasisStructure{Tensor} family. Need one vector per dimension
    if length(x) == N
        return x
    end

    # otherwise throw an error
    error("Basis is $N dimensional, need one Vector per dimension")
end

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

# code to be run at the top of each `BasisStructure` constructor
# it enforces compatibility of arguments and computes common items
function check_basis_structure{N}(basis::Basis{N}, x, order)
    if N > 1 && size(order, 2) == 1  # 62
        order = order * ones(Int, 1, N)  # 63
    end  # 64

    if N == 1 && isa(order, Int)
        order = fill(order, 1, 1)
    end

    # initialize basis structure (66-74)
    m = size(order, 1)
    if m > 1
        minorder = fill(minimum(order), 1, N)
        numbases = (maximum(order) - minorder) + 1
    else
        minorder = order + zeros(Int, 1, N)
        numbases = fill(1, 1, N)
    end

    x = _checkx(N, x)

    return m, order, minorder, numbases, x
end

# give the type of the `vals` field based on the family type parameter of the
# corresponding basis. `Spline` and `Lin` use sparse, `Cheb` uses dense
# a hybrid must fall back to a generic AbstractMatrix{Float64}
_vals_type(::Type{Spline}) = Base.SparseMatrix.SparseMatrixCSC{Float64,Int}
_vals_type(::Type{Lin}) = Base.SparseMatrix.SparseMatrixCSC{Float64,Int}
_vals_type(::Type{Cheb}) = Matrix{Float64}
_vals_type(::Type{BasisFamily}) = AbstractMatrix{Float64}

# conveneince method so we can pass an instance of the type also
_vals_type{TF<:BasisFamily}(::TF) = _vals_type(TF)

# --------------- #
# convert methods #
# --------------- #

# no-op
Base.convert{T<:ABSR}(::Type{T}, bs::BasisStructure{T}) = bs

# funbconv from direct to expanded
function Base.convert{TM}(bst::Type{Expanded}, bs::BasisStructure{Direct,TM},
                      order=fill(0, 1, size(bs.order, 2)))
    d, numbas, d1 = check_convert(bst, bs, order)
    n = prod([size(bs.vals[1, j], 2) for j=1:d])

    vals = Array(TM, numbas)

    for i=1:numbas
        vals[i] = bs.vals[order[i, d] - bs.order[d]+1, d]  # 63
        for j=d-1:-1:1
            vals[i] = row_kron(vals[i],
                               bs.vals[order[i, j] - bs.order[j]+1, j])  #65
        end
    end
    BasisStructure{Expanded,TM}(order, vals)
end

# funbconv from tensor to expanded
function Base.convert{TM}(bst::Type{Expanded}, bs::BasisStructure{Tensor,TM},
                      order=fill(0, 1, size(bs.order, 2)))
    d, numbas, d1 = check_convert(bst, bs, order)
    m = prod([size(bs.vals[1, j], 1) for j=1:d])
    n = prod([size(bs.vals[1, j], 2) for j=1:d])

    vals = Array(TM, numbas)

    for i=1:numbas  # 54
        vals[i] = bs.vals[order[i, d] - bs.order[d]+1, d]  # 55
        for j=d-1:-1:1  # 56
            # 57
            vals[i] = kron(vals[i], bs.vals[order[i, j] - bs.order[j]+1, j])
        end
    end

    BasisStructure{Expanded,TM}(order, vals)
end

# funbconv from tensor to direct
# TODO: there is probably a more efficient way to do this, but since I don't
#       plan on doing it much, this will do for now. The basic point is that
#       we need to expand the rows of each element of `vals` so that all of
#       them have prod([size(v, 1) for v in bs.vals])) rows.
function Base.convert{TM}(bst::Type{Direct}, bs::BasisStructure{Tensor,TM},
                      order=fill(0, 1, size(bs.order, 2)))
    d, numbas, d1 = check_convert(bst, bs, order)
    vals = Array(TM, numbas, d)
    raw_ind = cell(1, d)  # 78

    for j=1:d
        for i=1:size(bs.vals, 1)
            if !(isempty(bs.vals[i, j]))  # 84
                raw_ind[j] = collect(1:size(bs.vals[i, j], 1))  # 85
                break
            end
        end
    end

    ind = gridmake(raw_ind...)  # 90

    for j=1:d, i=1:numbas
        if !(isempty(bs.vals[i, j]))  # 93
            vals[i, j] = bs.vals[i, j][ind[:, j], :]  # 94
        end
    end

    BasisStructure{Direct,TM}(order, vals)
end

# ------------ #
# Constructors #
# ------------ #

# quick function to take order+vals and return expanded form for 1d problems
function to_expanded(out_order::Matrix{Int}, vals::Array)
    vals = vals[collect(out_order + (1 - minimum(out_order)))]
    BasisStructure{Expanded,eltype(vals)}(out_order, vals)
end


# method to construct BasisStructure in direct or expanded form based on
# a matrix of `x` values  -- funbasex
function BasisStructure{N,BF}(basis::Basis{N,BF}, ::Direct,
                              x::Array{Float64}=nodes(basis)[1], order=0)

    m, order, minorder, numbases, x = check_basis_structure(basis, x, order)
    # 76-77
    out_order = minorder
    out_format = Direct()
    val_type = _vals_type(BF)
    vals = Array(val_type, maximum(numbases), N)

    # now do direct form, will convert to expanded later if needed
    for j=1:N
        # 126-130
        if (m > 1)
            orderj = unique(order[:, j])
        else
            orderj = order[1, j]
        end

        #131-135
        if length(orderj) == 1
            vals[1, j] = evalbase(basis.params[j], x[:, j], orderj)[1]
        else
            vals[orderj-minorder[j]+1, j] =
                evalbase(basis.params[j], x[:, j], orderj)[1]
        end
    end

    # if N == 1, switch to expanded format and return it directly
    # 140-145
    if N == 1  # 1 dimension
        return to_expanded(order, vals)
    end

    # construct Direct Format
    BasisStructure{Direct,val_type}(out_order, vals)
end

function BasisStructure(basis::Basis, ::Expanded,
                        x::Array{Float64}=nodes(basis)[1], order=0)  # funbasex
    # create direct form, then convert to expanded
    bsd = BasisStructure(basis, Direct(), x, order)
    convert(Expanded, bsd)
end

function BasisStructure{N,BF,T}(basis::Basis{N,BF}, ::Tensor,
                           x::Vector{Vector{T}}=nodes(basis)[2], order=0)
    m, order, minorder, numbases, x = check_basis_structure(basis, x, order)
    out_order = minorder
    out_format = Tensor()
    val_type = _vals_type(BF)
    vals = Array(val_type, maximum(numbases), N)

    # construct tensor base
    for j=1:N
        # 113-117
        if (m > 1)
            orderj = unique(order[:, j])
        else
            orderj = order[1, j]
        end

        #118-122
        if length(orderj) == 1
            vals[1, j] = evalbase(basis.params[j], x[j], orderj)[1]
        else
            vals[orderj-minorder[j]+1, j] = evalbase(basis[:params][j], x[j],
                                                     orderj)[1]
        end
    end

    N == 1 ? to_expanded(order, vals) :
             BasisStructure{Tensor,val_type}(out_order, vals)
end
