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

    @assert nobsa == nobsb "A and B must have same number of rows"

    # fill in each element. To do this we make sure we access each array
    # consistent with its column major memory layout.
    for ia=1:na, ib=1:nb, t=1:nobsa
        @inbounds out[t, nb*(ia-1) + ib] = A[t, ia] * B[t, ib]
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
    out = spzeros(promote_type(S, T), nobsa, na*nb)
    row_kron!(A, B, out)
    out
end

const dprod = row_kron

# cckronxi.m -- DONE
cckronxi(b::Matrix, c, ind=1:length(b)) = b \ c  # 23

function cckronxi(b::Vector{Any}, c, ind=1:length(b))
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

# cdprodx.m -- DONE
cdprodx{T<:Number}(b::Matrix{T}, c, ind=1:prod(size(b))) = b*c  # 39


# TODO: this should be a fold
function cdprodx(b::Array{Any}, c, ind=1:prod(size(b)))
    d = length(ind)
    a = b[ind[d]]
    for i=d-1:-1:1
        a = dprod(b[ind[i]], a)
    end
    a = a * c
end


# cckronx.m -- DONE
cckronx(b::Matrix, c, ind=1:prod(size(b))) = b * c  # 23

function cckronx(b::Vector{Any}, c, ind=1:prod(size(b)))
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
