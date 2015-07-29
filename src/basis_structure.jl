# ------------------- #
# BasisStructure Type #
# ------------------- #

abstract AbstractBasisStructureRep
typealias ABSR AbstractBasisStructureRep

immutable Tensor <: ABSR end
immutable Direct <: ABSR end
immutable Expanded <: ABSR end

immutable BasisStructure{BST<:ABSR}
    order::Matrix{Int}
    vals::Array{AbstractMatrix}
end

Base.writemime{BST}(io::IO, ::MIME"text/plain", b::BasisStructure{BST}) =
    print("BasisStructure{$BST} of order $(b.order)")

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

# no-op
Base.convert{T<:ABSR}(::Type{T}, bs::BasisStructure{T}) = bs

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
                      order=fill(0, 1, size(bs.order, 2)))
    d, numbas, d1 = check_convert(bst, bs, order)
    vals = Array(AbstractMatrix, numbas, d)
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

    BasisStructure{Direct}(order, vals)
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

# TODO: this method is clobbered by the one for bformat=Tensor() below.
#       need to decide if I want default to be Direct or Expanded. It would
#       be much easier to do this if bformat were the second argument so
#       that I could have `BasisStucture(basis)` give default and
#       `BasisStructure(basis, bformat)` give alternate
# method to construct BasisStructure in direct or expanded form based on
# a matrix of `x` values  -- funbasex
function BasisStructure{N}(basis::Basis{N}, ::Direct,
                           x::Array{Float64}=nodes(basis)[1], order=0)

    d, m, order, minorder, numbases = check_basis_structure(basis, x, order,
                                                            Direct())
    # 76-77
    out_order = minorder
    out_format = Direct()
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
            vals[1, j] = evalbase(basis.params[j], x[:, j], orderj)[1]
        else
            vals[orderj-minorder[j]+1, j] =
                evalbase(basis.params[j], x[:, j], orderj)[1]
        end
    end

    # if d == 1, switch to expanded format and return it directly
    # 140-145
    if N == 1  # 1 dimension
        return to_expanded(order, vals)
    end

    # construct Direct Format
    BasisStructure{Direct}(out_order, vals)
end

function BasisStructure(basis::Basis, ::Expanded,
                        x::Array{Float64}=nodes(basis)[1], order=0)  # funbasex
    # create direct form, then convert to expanded
    bsd = BasisStructure(basis, Direct(), x, order)
    convert(Expanded, bsd)
end

function BasisStructure{N,T}(basis::Basis{N}, ::Tensor,
                           x::Vector{Vector{T}}=nodes(basis)[2], order=0)
    d, m, order, minorder, numbases = check_basis_structure(basis, x, order,
                                                            Tensor())
    out_order = minorder
    out_format = Tensor()
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
            vals[1, j] = evalbase(basis.params[j], x[j], orderj)[1]
        else
            vals[orderj-minorder[j]+1, j] = evalbase(basis[:params][j], x[j],
                                                     orderj)[1]
        end
    end

    N == 1 ? to_expanded(order, vals) : BasisStructure{Tensor}(out_order, vals)
end

# add method to funbasex that calls the above constructors
funbasex(basis::Basis, x=nodes(basis)[1], order=0, bformat::ABSR=Direct()) =
    BasisStructure(basis, x, order, bformat)

funbase(basis::Basis, x=nodes(basis)[1], order=fill(0, 1, ndims(basis))) =
    funbasex(basis, x, order, Expanded()).vals[1]
