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

fix{T<:Real}(x::Array{T}) = fix!(x, similar(x, Int))
fix{T<:Real}(x::T) = x >= 0 ? floor(Int, x) : ceil(Int, x)

# ckronx.m -- DONE
function ckronx{TM<:AbstractMatrix}(b::Matrix{TM}, c::Array,
                                    ind::Matrix{Int}=reshape(1:length(b),
                                                             1, length(b)))
    d = length(ind)  # 26
    n = Array(Int, d)  # 27
    for i=1:d  # 28
        n[i] = size(b[ind[i]], 2)
    end

    if prod(n) != size(c, 1)  # 29-31
        m = "b and c are not conformable (b suggests size(c, 1) should be $(prod(n)))"
        error(m)
    end

    z = c'  # 32
    mm = 1  # 33
    for i=1:d
        @compat m = Int(length(z) / n[i])  # 35
        z = reshape(z, m, n[i])  # 36
        z = b[ind[i]] * z'  # 37
        mm = mm * size(z, 1)  # 38
    end
    z = reshape(z, mm, size(c, 2))  # 40
end

# ckron.m -- DONE
ckron(A::AbstractArray, B::AbstractArray) = kron(A, B)
ckron(arrays::AbstractArray...) = reduce(kron, arrays)

gridmake(v::AbstractVector) = v

# gridmake.m -- DONE
function gridmake{T}(arrays::AbstractVector{T}...)
    # TODO: this gridmake works, but I don't like it.
    shapes = Int[size(e, 1) for e in arrays]

    n = length(arrays)
    l = prod(shapes)
    out = Array(T, l, n)

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
        a = row_kron(b[ind[i]], a)
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

# lookup.c -- DONE
function lookup(table::Vector, x::Vector, p::Int=0)
    n = length(table)
    m = length(x)
    out = fill(42, m)

    # lower enbound adjustment
    numfirst = 1
    t1 = table[1]
    for i=2:n
        if table[i] == t1
            numfirst += 1
        else
            break
        end
    end

    # upper endpoint adjustment
    tn = table[n]
    if p >= 2
        n -= 1
        for i=n:-1:1
            if table[i] == tn
                n -= 1
            else
                break
            end
        end
    end

    n1 = n - 1
    n2 = n - 1

    if n - numfirst < 1  # only one unique value in table
        if p == 1 || p == 3
            for i=1:m
                out[i] = numfirst
            end
        else
            for i=1:m
                if table[1] <= x[i]
                    out[i] = numfirst
                else
                    out[i] = 0
                end
            end
        end
        return out
    end

    jlo = 1

    for i=1:m
        inc = 1
        xi = x[i]
        if xi >= table[jlo]
            jhi = jlo + 1
            while jhi <= n && xi >= table[jhi]
                jlo = jhi
                jhi += inc
                if jhi > n
                    jhi = n+1
                end
            end
        else
            jhi = jlo
            jlo -= 1
            while jlo > 0 && xi < table[jlo]
                jhi = jlo
                jlo -= inc
                if jlo < 1
                    jlo = 0
                    break
                else
                    inc += inc
                end
            end
        end

        while jhi - jlo > 1
            j = (jhi + jlo) >> 1
            if xi >= table[j]
                jlo = j
            else
                jhi = j
            end
        end

        out[i] = jlo

        if jlo < 1
            jlo = 1
            if p == 1 || p == 3
                out[i] = numfirst
            end
        end

        if jlo == n1
            jlo = n2
        end
    end

    out

end


# utility to expand the order input if needed
# used in basis_structure.jl and interp.jl
_check_order(N::Int, order::Int) = fill(order, 1, N)
_check_order(N::Int, order::Vector) = reshape(order, 1, N)

function _check_order(N::Int, order::Matrix)
    if size(order, 2) == N
        return order
    end

    if size(order, 1) == N
        m = size(order, 2)
        return reshape(order, m, N)
    end

    error("invalid order argument. Should have $N columns")
end
