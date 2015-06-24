# ---------------- #
# Helper functions #
# ---------------- #

# fix.m -- DONE
function fix!{T <: Real}(x::Array{T}, out::Array{Int})
    for i=1:length(x)  # use linear indexing
        out[i] = fix(x[i])
    end
    return out
end
fix{T <: Real}(x::Array{T}) = fix!(x, similar(x, Int))
fix{T <: Real}(x::T) = int(x >= 0 ? floor(x) : ceil(x))

# ckron.m -- DONE
ckron(A::Array, B::Array) = kron(A, B)
ckron(arrays::Array...) = reduce(kron, arrays)

gridmake(v::Vector) = v

# gridmake.m -- DONE
function gridmake(arrays::Vector...)
    # TODO: this gridmake works, but I don't like it.
    shapes = Int[size(e, 1) for e in arrays]

    n = length(arrays)
    l = prod(shapes)
    out = Array(Float64, l, n)

    shapes = shapes[end:-1:1]
    sh = push!([1], shapes[1:end-1]...)
    repititions = cumprod(sh)
    repititions = repititions[end:-1:1]

    for i=1:n
        arr = arrays[i]
        outer = repititions[i]
        inner = floor(Int, l / (outer * size(arr, 1)))
        out[:, i] = repeat(arrays[i], inner=[inner], outer=[outer])
    end
    return out
end

# dprod.m  - DONE
function row_kron!(A::AbstractMatrix, B::AbstractMatrix, out::AbstractMatrix)
    # get input dimensions
    nobsa, na = size(A)
    nobsb, nb = size(B)

    nobsa != nobsb && error("A and B must have same number of rows")

    # fill in each element. To do this we make sure we access each array
    # consistent with its column major memory layout.
    @inbounds for ia=1:na, ib=1:nb, t=1:nobsa
        out[t, nb*(ia-1) + ib] = A[t, ia] * B[t, ib]
    end
    out
end

function row_kron{S,T}(A::AbstractMatrix{S}, B::AbstractMatrix{T})
    nobsa, na = size(A)
    nobsb, nb = size(B)
    out = Array(promote_type(S, T), nobsa, na*nb)
    row_kron!(A, B, out)
    out
end

function row_kron{S,T}(A::SparseMatrixCSC{S}, B::SparseMatrixCSC{T})
    nobsa, na = size(A)
    nobsb, nb = size(B)

    # doing this on the transpose so the row indices will be sorted
    cols_a, rows_a, vals_a = findnz(A')
    cols_b, rows_b, vals_b = findnz(B')

    # nnza, nnzb = map(length, (ra, rb))

    prev_last_a = searchsortedfirst(rows_a, 0)
    prev_last_b = searchsortedfirst(rows_b, 0)

    I = Array(Int64, 0)
    J = Array(Int64, 0)
    V = Array(promote_type(S, T), 0)

    for t in 1:nobsa
        next_last_a = searchsortedfirst(rows_a, t+1)
        next_last_b = searchsortedfirst(rows_b, t+1)
        these_cols_a = cols_a[prev_last_a:next_last_a-1]
        these_cols_b = cols_b[prev_last_b:next_last_b-1]

        these_vals_a = vals_a[prev_last_a:next_last_a-1]
        these_vals_b = vals_b[prev_last_b:next_last_b-1]


        for ia in 1:length(these_cols_a)
            ca = these_cols_a[ia]
            for ib in 1:length(these_cols_b)
                cb = these_cols_b[ib]
                push!(I, t)
                push!(J, nb*(ca-1) + cb)
                push!(V, these_vals_a[ia] * these_vals_b[ib])
            end
        end


        prev_last_a = next_last_a
        prev_last_b = next_last_b

        # # cut down ra, rb so search sorted doesn't have to work as hard
        # rows_a = rows_a[prev_last_a-1:end]
        # rows_b = rows_b[prev_last_a-1:end]
    end

    sparse(I, J, V, nobsa, na*nb)
end

const dprod = row_kron

# ckronxi.m -- DONE
ckronxi{T<:Number}(b::Matrix{T}, c, ind=1:length(b)) = b \ c  # 23

function ckronxi(b::Array, c, ind=1:length(b))
    d = length(ind)  # 25
    n = Int[size(b[ind[i]], 2) for i=1:d]  #26-27
    prod(n) != size(c, 1) && error("b and c are not conformable")  # 28-30

    z = c'  # 31
    mm = 1  # 32
    for i=1:d  # 33
        m = round(Int, prod(size(z)) / n[i])  # 34
        z = reshape(z, m, n[i])  # 35
        z = b[ind[i]] \ z'  # 36
        mm *= size(z, 1)  # 37
    end  # 38
    reshape(z, mm, size(c, 2))  # 39
end

# cdprodx.m -- DONE
cdprodx{T<:Number}(b::Matrix{T}, c, ind=1:prod(size(b))) = b*c  # 39


# TODO: this should be a fold
function cdprodx(b::Array, c, ind=1:prod(size(b)))
    d = length(ind)
    a = b[ind[d]]
    for i=d-1:-1:1
        a = dprod(b[ind[i]], a)
    end
    a = a * c
end


# cckronx.m -- DONE
cckronx{T<:Number}(b::Matrix{T}, c, ind=1:prod(size(b))) = b * c  # 23

function cckronx(b::Array, c, ind=1:prod(size(b)))
    d = length(ind)  # 25
    n = Int[size(b[ind[i]], 2) for i=1:d]  #26-27
    prod(n) != size(c, 1) && error("b and c are not conformable")  # 28-30

    z = c'  # 31
    mm = 1  # 32
    for i=1:d  # 33
        m = prod(size(z)) / n[i]  # 34
        z = reshape(z, m, n[i])  # 35
        z = b[ind[i]] \ z'  # 36
        mm = mm*size(z, 1)  # 37
    end  # 38
    reshape(z, mm, size(c, 2))  # 39
end

# nodeunif.m -- DONE
function nodeunif(n::Int, a::Int, b::Int)
    x = linspace(a, b, n)
    return x, x
end

function nodeunif(n::Array, a::Array, b::Array)
    d = length(n)
    xcoord = cell(d)
    for k=1:d
        xcoord[k] = linspace(a[k], b[k], n[k])
    end
    return gridmake(xcoord...), xcoord
end

function squeeze_trail(x::Array)
    sz = size(x)
    squeezers = Int[]
    n = length(sz)
    for i=n:-1:1
        if sz[i] == 1
            push!(squeezers, i)
        else
            break
        end
    end
    squeeze(x, tuple(squeezers...))
end


# lookup.m -- DONE
function lookup(tabvals::Vector, x::Vector, endadj=0)
    n = prod(size(x))
    m = length(tabvals)
    if endadj >= 2
        m = m - sum(tabvals .== tabvals[end])
    end

    ind = sortperm(vcat(tabvals[1:m], x))
    temp = find(ind .>m)
    j = ind[temp] - m
    ind = reshape(temp .- (1:n), size(x)...)
    ind[j] = ind[:]

    if endadj == 1 || endadj == 3
        ind[ind .== 0] = sum(tabvals .== tabvals[1])
    end
    ind
end
