# TODO: still need to write fund, minterp
# TODO: also need splidop, lindop
# TODO: funeval fails for scalar input and does weird thing for 1-element
#       vector input


const BASE_TYPES = [:spli, :cheb, :lin]

base_exists(s::Symbol) = s in BASE_TYPES

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




# -------------------- #
# Dispatcher functions #
# -------------------- #


basedef(s::Symbol, args...) =
    s == :spli ? splidef(args...) :
    s == :cheb ? chebdef(args...) :
    s == :lin  ? lindef(args...)  :
    error("somehow you snuck through here you :japanese_goblin:")

basenode(s::Symbol, args...) =
    s == :spli ? splinode(args...) :
    s == :cheb ? chebnode(args...) :
    s == :lin  ? linnode(args...)  :
    error("somehow you snuck through here you :japanese_goblin:")

evalbase(s::Symbol, args...) =
    s == :spli ? splibase(args...) :
    s == :cheb ? chebbase(args...) :
    s == :lin  ? linbase(args...)  :
    error("somehow you snuck through here you :japanese_goblin:")


# ---------------------------- #
# Generic translated functions #
# ---------------------------- #

# from fundef.m  -- DONE
function fundef(foo...)
    d = length(foo)  # 89
    n = zeros(Int, d)  # 93
    b = zeros(d)  # 94
    a = zeros(d)  # 95
    p = cell(d)  # 96

    basetype = Array(Symbol, d)
    for j=1:d
        basetype[j] = foo[j][1]  # 99
        !(base_exists(basetype[j])) && error("Unknown basis $(foo[j][1])")
        n[j], a[j], b[j], p[j] = basedef(basetype[j], foo[j][2:end]...)  # 124
    end

    # package output. Lines 143-150
    g = Dict{Symbol, Any}()
    g[:d] = d
    g[:n] = n
    g[:a] = a
    g[:b] = b
    g[:basetype] = basetype
    g[:parms] = p
    g
end


# fundefn.m
function fundefn(basistype::Symbol, n, a, b, order=3)
    d = length(n)
    length(a) != d && error("a must be same dimension as n")
    length(b) != d && error("b must be same dimension as n")
    any(a .> b) && error("left endpoints must be less than right endpoints")
    any(n .< 2) && error("n(i) must be greater than 1")

    params = cell(1, d)
    if basistype == :cheb
        for i=1:d params[i] = Any[:cheb, n[i], a[i], b[i]] end
    elseif basistype == :spli
        for i=1:d params[i] = Any[:spli, [a[i], b[i]], n[i]-order+1, order] end
    elseif basistype == :lin
        for i=1:d params[i] = Any[:lin, [a[i], b[i]], n[i]] end
    end

    fundef(params...)
end

# funnode.m -- DONE
function funnode(basis::Dict)
    d = basis[:d]  # 18
    if d == 1  # 19
        xcoord = basenode(basis[:basetype][1], basis[:parms][1]...)  # 20
        x = xcoord  # 21
    else  # 22
        xcoord = cell(d)  # 23
        for j=1:d  # 24
            xcoord[j] = basenode(basis[:basetype][j], basis[:parms][j]...)  #25
        end  # 26
        x = gridmake(xcoord...)  # 27
    end

    return x, xcoord
end

# funbase.m -- DONE
funbase(basis::Dict, x=funnode(basis)[1],
        order=fill(0, 1, basis[:d])) =
    funbasex(basis, x, order, :expanded)[:vals][1]


# funbasex.m -- DONE
function funbasex(basis, x=funnode(basis)[1], order=0, bformat::Symbol=:none)
    d = length(basis[:n])  # 59
    if d > 1 && size(order, 2) == 1  # 62
        order = order * ones(Int, 1, d)  # 63
    end  # 64

    # initialize basis structure (66-74)
    m = size(order, 1)
    if m > 1
        minorder = fill(minimum(order), 1, d)
        numbases = (maximum(order) - minorder) + 1
    else
        minorder = order + zeros(Int, 1, d)
        numbases = fill(1, 1, d)
    end

    B = Dict{Symbol, Any}()  # 76
    B[:order] = minorder
    B[:format] = bformat
    B[:vals] = cell(maximum(numbases), d)  # 77

    # 83-89
    if bformat == :none
        if isa(x, Array{Any})
            bformat = :tensor
        else
            bformat = :direct
        end
    end

    # 91-101
    if d > 1
        if !(isa(x, Array{Any})) && bformat == :tensor
            error("Must pass a cell array to form a tensor format basis")
        end
        if isa(x, Array{Any}) && bformat == :direct
            # it would be more efficient in this case to
            # use the tensor form to compute the bases and then
            # to use indexing to expand to the direct form
            #      x=gridmake(x); % convert to grid for direct form
        end
    end

    #103-107
    B[:format] = is(x, Array{Any}) ? :tensor : :direct

    if B[:format] == :tensor  # 111
        for j=1:d
            # 113-117
            if (m > 1)
                orderj = unique(order[:, j])
            else
                orderj = order[1, j]
            end

            #118-122
            if length(orderj) == 1
                B[:vals][1, j] = evalbase(basis[:basetype][j],
                                          basis[:parms][j]..., x[j], orderj)
            else
                B[:vals][orderj-minorder[j]+1, j] =
                    evalbase(basis[:basetype][j], basis[:parms][j]..., x[j],
                             orderj)
            end
        end
    else  # B[:format] == :direct
        for j=1:d
            # 126-130
            if (m > 1)
                orderj = unique(order[:, j])
            else
                orderj = order[1, j]
            end

            #131-135
            if length(orderj) == 1
                B[:vals][1, j] = evalbase(basis[:basetype][j],
                                          basis[:parms][j]..., x[:, j], orderj)[1]
            else
                B[:vals][orderj-minorder[j]+1, j] =
                    evalbase(basis[:basetype][j], basis[:parms][j]..., x[:, j],
                             orderj)[1]
            end
        end
    end

    # if d == 1, switch to expanded format
    # 140-145
    if size(B[:vals], 2) == 1  # 1 dimension
        B[:format] = :expanded
        B[:order] = order
        B[:vals] = B[:vals][collect(order + (1 - minimum(order)))]
        return B
    end

    # create expanded format  148-153
    if bformat == :expanded
        B = funbconv(B, order, :expanded)
    elseif bformat == :direct
        if isa(x, Array{Any})
            B = funbconv(B, order, :direct)
        end
    end
    return B
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

# funfitf.m -- DONE
function funfitf(basis, f::Function, args...)
    x = funnode(basis)[1]
    y = f(x, args...)
    funfitxy(basis, x, y)[1]
end

# funfitxy.m -- DONE
function funfitxy(basis, x, y)
    m = size(y, 1)
    prod(basis[:n]) > m && error("Cannot be more basis functions than points")

    if isa(x, Dict)  # 61 : use precomputed basis structure
        B = x

        if B[:format] == :tensor
            if any(B[:order][1, :] != 0)
                m = "invalid basis structure - first elements must be order 0"
                error(m)
            end
            B[:vals] = B[:vals][1, :]  # 68
            c = ckronxi(B[:vals], y, basis[:d]:-1:1)  # 66

        elseif B[:format] == :direct
            B = funbconv(B, zeros(1, size(B[:vals], 2)), :expanded)  # 71
            c = B[:vals][1] \ y  # 72

        elseif B[:format] == :expanded
            c = B[:vals][1] \ y  # 74
        end

    elseif isa(x, Array{Any})  # evaluate at grid points
        # 77
        mm = 1
        for i=1:size(x, 2)
            mm = mm * size(x[i], 1)
        end

        mm != m && error("x and y are incompatible")

        B = funbasex(basis, x, 0)  # 81
        c = ckronxi(B[:vals], y, basis[:d]:-1:1)  # 82

    else  # evaluate at arbitrary points
        size(x, 1) != m && error("x and y are incompatible")

        B = funbasex(basis, x, 0, :expanded)
        c = B[:vals][1] \ y
    end
    return c, B
end


# funeval.m -- DONE
function funeval(c, basis, B, order=0)
    isempty(c) && error("missing basis coefficients")

    if !isa(B, Dict)  # B is not a basis structure, so construct one
        d = basis[:d]  # 60
        if size(B, 2) != d  # TODO: this will be ill-defiend for non array B's
            error("x must have d=$(d) columns")  # 62
        end

        if size(order, 2) == 1
            order = order * fill(1, 1, d)  # 65
        end
        B = funbasex(basis, B, order)  # 67
    else
        if size(order, 2) == 1
            order = order * ones(Int, 1, size(B[:order], 2))  # 69
        end
    end

    # now we we have a basis structure
    if isa(B[:vals], Array{Any})
        y = B[:format] == :tensor   ? funeval1(c, B, order) :  # 75
            B[:format] == :direct   ? funeval2(c, B, order) :  # 77
            B[:format] == :expanded ? funeval3(c, B, order) :  # 79
            error("Basis structure has invalid ``format`` field")  # 81
    else
        y = B[:vals] * c  # 84
    end

    return y, B
end

# inside funeval.m. Evaluates according to ``:tensor`` format -- DONE
function funeval1(c, B, order=fill(0, 1, size(B[:order], 2)))
    kk, d = size(order)  # 95
    # 98 reverse the order of evaluation: B(d)xB(d-1)x...xB(1)
    order = flipdim(order+1*(size(B[:vals], 1) * (0:d-1)' - B[:order]+1), 2)

    # 99
    nx = 1
    for j=1:d
        nx = nx * size(B[:vals][i, j], 1)
    end

    f = zeros(nx, size(c, 2), kk)  # 100

    for i=1:kk
        f[:, :, i] = ckronx(B.vals, c, order[i, :])  # 102
    end
    return squeeze_trail(f)
end

# inside funeval.m. Evaluates according to ``:direct`` format -- DONE
function funeval2(c, B, order=fill(0, 1, size(B[:order], 2)))
    kk, d = size(order)  # 95
    # 114 reverse the order of evaluation: B(d)xB(d-1)x...xB(1)
    order = flipdim(order+1*(size(B[:vals], 1) * (0:d-1)' - B[:order]+1), 2)

    f = zeros(size(B[:vals][1], 1), size(c, 2), kk)  # 116

    for i=1:kk
        f[:, :, i] = cdprodx(B[:vals], c, order[i, :])  # 118
    end
    return squeeze_trail(f)
end

# inside funeval.m. Evaluates according to ``:expanded`` format -- DONE
function funeval3(c, B, order=fill(0, 0, 0))
    if isempty(order)
        if isa(B[:vals], Array{Any})
            kk = length(B[:vals])  # 126
            order = 1:kk'  # 127
        else
            kk = 1  # 128
            order = 1  # 129
        end
    else
        kk = size(order, 1)  # 133
    end

    if isa(B[:vals], Array{Any})
        nx = size(B[:vals][1], 1)
        f = zeros(nx, size(c, 2), kk)
        for i=1:kk
            # 140 determine which element of B[:vals] is the desired basis
            ii = Int[]
            for row=1:size(B[:order], 1)
                r = B[:order][row, :]
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

            f[:, :, i] = B[:vals][ii]*c  #148
         end
     else
        nx = size(B[:vals], 1)  # 151
        f = zeros(nx, size(c, 2), kk)  # 152
        for i=1:kk
            f[:, :, i] = B[:vals]*c  # 154
        end
    end

    return squeeze_trail(f)
end

# fund.m
function fund(c, basis, x, hess_opt)
    # TODO: come back when I need this. I think I should probably do something
    #       like what optim does and write functions `f`, `fg!` and `fgh!` to
    #       replicate this for the type of basis instead of this function
    nothing
end




# funbconv.m  -- DONE
function funbconv(b, order=fill(0, 1, size(b[:order], 2)),
                  format::Symbol=:expanded)

    d = size(b[:order], 2)  # 24
    numbas, d1 = size(order, 1), size(order, 2)  # 34
    d1 != d && error("ORDER incompatible with basis functions")  # 35-37

    # 39-41
    if any(minimum(order) .< b[:order])
        error("Order of derivative operators exceeds basis")
    end

    B = Dict{Symbol, Any}()
    if format == :expanded
        B[:vals] = cell(numbas)
        B[:order] = order
        B[:format] = format

        if b[:format] == :tensor
            m, n = 1, 1
            for j=1:d  # 50
                m = m * size(b[:vals][1, j], 1)  # 51
                n = n * size(b[:vals][1, j], 2)  # 52
            end

            for i=1:numbas  # 54
                B[:vals][i] = b[:vals][order[i, d] - b[:order][d]+1, d]  # 55
                for j=d-1:-1:1  # 56
                    # 57
                    B[:vals][i] = kron(B[:vals][i], b[:vals][order[i, j] - b[:order][j]+1, j])
                end
            end
        elseif b[:format] == :direct
            n = 1  # 61
            for j=1:d  # 61
                n = n * size(b[:vals][1, j], 2)  # 61
            end

            for i=1:numbas
                B[:vals][i] = b[:vals][order[i, d] - b[:order][d]+1, d]  # 63
                for j=d-1:-1:1
                    B[:vals][i] = row_kron(B[:vals][i],
                                           b[:vals][order[i, j] - b[:order][j]+1, j])  #65
                end
            end
        elseif b[:format] == :expanded
            warn("Basis already in expanded form")
            B = b
        else
            error("Improper basis format")
        end

    elseif format == :direct && b[:format] == :tensor
        B[:vals] = cell(numbas, d)  # 75
        B[:order] = order  # 76
        B[:format] = format  # 77
        ind = cell(1, d)  # 78

        for j=1:d, i=1:size(b[:vals], 1)
            if !(isempty(b[:vals][i, j]))  # 81
                ind[j] = collect(1:size(b[:vals][i, j], 1))  # 82
                break
            end
        end

        ind = gridmake(ind...)  # 87

        for j=1:d, i=1:numbas
            if !(isempty(b[:vals][i, j]))  # 90
                B[:vals][i, j] = B[:vals][i, j][ind[j], :]
            end
        end
    else
        error("Not implemented for this option")
    end

    return B
end


# --------------------- #
# Cheb-family functions #
# --------------------- #

# chebdef.m -- DONE
function chebdef(n::Int, a::Real, b::Real)
    n <= 0 && error("n must be positive")
    a >= b && error("left endpoint (a) must be less than right end point (b)")
    return n, a, b, Any[n, a, b]
end


# chebnode.m -- DONE
function chebnode(n::Int, a::Real, b::Real, nodetype=0)
    s = (b-a) / 2  # 21
    m = (b+a) / 2  # 22

    if nodetype < 2  # usual nodes
        k = π*(0.5:(maximum(n) - 0.5))'  # 25
        x = m - cos(k[1:n] / n) * s  # 26
        if nodetype == 1  # Extend nodes to endpoints
            aa = x[1]  # 28
            bb = x[end]  # 29
            x = (bb*a - aa*b)/(bb-aa) + (b-a)/(bb-aa)*x  # 30
        end
    else # Lobatto nodes
        k = pi*(0:(maximum(n)-1))  # 33
        x = m - cos(k[1:n] / (n-1)) * s  # 34
    end
    x
end


# chebdop.m -- DONE
function chebdop(n, a, b, order=1)
    if order > 0
        # TODO: figure out some caching mechanism that avoids globals
        D = cell(max(2, order), 1)  # 49
        i = repeat(collect(1:n), outer=[1, n]); j = i'  # 50

        # 51
        inds = find((rem(i + j, 2) .== 1) & (j.>i))
        r, c = similar(inds), similar(inds)
        for i in 1:length(inds)
            r[i], c[i] = ind2sub((n, n), inds[i])
        end

        d = sparse(r, c, (4/(b-a)) * (vec(j[1, c])-1), n-1, n)  # 52
        d[1, :] ./= 2  # 53
        D[1] = d  # 54
        for ii=2:max(2, order)
            D[ii] = d[1:n-ii, 1:n-ii+1] * D[ii-1]  # 56
        end
    else
        D = cell(abs(order), 1)  # 64
        nn = n - order  # 65
        i = (0.25 * (b - a)) ./(1:nn)  # 66
        d = sparse(vcat(1:nn, 1:nn-2), vcat(1:nn, 3:nn), vcat(i, -i[1:nn-2]),
                   nn, nn)  # 67
        d[1, 1] *= 2  # 68
        d0 = ((-1).^(0:nn-1)') .* sum(d, 1)  # 69
        D[1] = sparse(vcat(d0[1:n]', d[1:n, 1:n]))  # 70
        for ii=-2:-1:order
            ind = 1:n-ii-1
            D[-ii] = sparse(vcat(d0[ind]', d[ind, ind]) * D[-ii-1])
        end
    end
    D, n-order, a, b, Any[n, a, b]
end


# chebbase.m -- DONE
function chebbase(n, a, b, x=chebnode(n, a, b, 1), order=0, nodetype=1)
    minorder = min(0, minimum(order))  # 30

    # compute 0-order basis
    if nodetype == 0
        temp = ((n-0.5):-1:0.5)''  # 41
        bas = cos((pi/n)*temp.*(0:(n-1-minorder))')  # 42
    else
        bas = chebbasex(n-minorder, a, b, x)  # 44
    end

    if length(order) == 1
        if order != 0
            D = chebdop(n, a, b, order)[1]
            B = bas[:, 1:n-order]*D[abs(order)]
        else
            B = bas
        end
    else
        B = cell(length(order))
        maxorder = maximum(order)
        if maxorder > 0 D = chebdop(n, a, b, maxorder)[1] end
        if maxorder < 0 I = chebdop(n, a, b, minorder)[1] end
        for ii=1:length(order)
            if order[ii] == 0
                B[ii] = bas[:, 1:n]
            elseif order[ii] > 0
                B[ii] = bas[:, i1:n-order[ii]] * D[order[ii]]
            else
                B[ii] = bas[:, i1:n-order[ii]] * I[-order[ii]]
            end
        end
    end

    return B, x
end

# chebbasex.m -- DONE
function chebbasex(n, a, b, x)
    z = (2/(b-a)) * (x-(a+b)/2)
    m = size(z, 1)
    bas = Array(Float64, m, n);
    bas[:, 1] = 1.0
    bas[:, 2] = z
    z = 2 * z
    for i=3:n
      bas[:, i] = z .* bas[:, i-1] - bas[:, i-2]
    end
    bas
end

# --------------------- #
# Spli-family functions #
# --------------------- #

# splidef.m -- DONE
function splidef(breaks, evennum=0, k::Int=3)
    # error handling
    k < 0 && error("spline order must be positive")
    length(breaks) < 2 && error("Must have at least two breakpoints")
    any(diff(breaks) .< 0) && error("Breakpoints must be non-decreasing")

    if evennum == 0  # 43
        if length(breaks) == 2  # 44
            evennum = 2
        end
    else
        if length(breaks) == 2
            breaks = linspace(breaks[1], breaks[2], evennum)
        else
            error("Breakpoint squence must contain 2 values when evennum > 0")
        end
    end

    n = length(breaks) + k - 1
    a = breaks[1]
    b = breaks[end]
    return n, a, b, Any[breaks, evennum, k]
end

# splinode.m  -- DONE
function splinode(breaks::Vector, evennum::Int, k::Int=3)
    a = breaks[1]  # 20
    b = breaks[end]  # 21
    n = length(breaks) + k - 1  # 22
    x = cumsum(vcat(fill(a, k), breaks, fill(b, k)))  # 23
    x = (x[1+k:n+k] - x[1:n]) / k  # 24
    x[1] = a  # 25
    x[end] = b  # 26
    x
end


# splidop.m
function splidop(breaks, evennum=0, k=3, order=1)
    any(order .>= k) && error("Order of differentiation must be less than k")

    # 38-40
    n = length(breaks) + k - 1
    kk = max(k - 1, k - order - 1)
    augbreaks = vcat(fill(breaks[1], kk), breaks, fill(breaks[end], kk))

    D = cell(abs(order), 1)
    if order > 0  # derivative
        temp = k ./ (augbreaks[k+1:n+k-1] - augbreaks[1:n-1])
        D[1] = spdiags([-temp temp], 0:1, n-1, n)  # TODO: pick up here
    end
end


# splibas.m -- DONE
function splibase(breaks::Vector, evennum, k=3, x=splinode(breaks, evennum, k),
                  order=0)
    # error handling
    k < 0 && error("spline order must be positive")
    !(issorted(breaks)) && error("Breakpoints must be non-decreasing")
    any(order .>= k) && error("Order of differentiation must be less than k")
    size(x, 2) > 1 && error("x must be a column vector")

    p = length(breaks)  # 53
    m = size(x, 1)  # 54
    minorder = minimum(order)  # 55

    # Augment the breakpoint sequence 57-59
    n = length(breaks)+k-1
    a = breaks[1]
    b = breaks[end]
    augbreaks = vcat(fill(a, k-minorder), breaks, fill(b, k-minorder))

    ind = lookup(augbreaks, x, 3)  # 69

    bas = zeros(m, k-minorder+1)  # 73
    bas[:, 1] = 1.0  # 74
    B = cell(length(order))  # 75

    # 76
    if maximum(order) > 0
        D = splidop(breaks, evennum, k, maximum(order))
    end

    # 77
    if minorder < 0
        I = splidop(breaks, evennum, k, minorder)
    end

    for j=1:k-minorder  # 78
        for jj=j:-1:1  # 79
            b0 = augbreaks[ind+jj-j]  # 80
            b1 = augbreaks[ind+jj]  # 81
            temp = bas[:, jj] ./ (b1 - b0)  # 82
            bas[:, jj+1] = (x - b0) .* temp + bas[:, jj+1]  # 83
            bas[:, jj] = (b1-x) .* temp  # 84
        end
        # bas not contains the order `j` spline basis
        ii = find((k-j) .== order)  # 87
        if !(isempty(ii))  # 88
            ii = ii[1]  # 89

            # Put values in appropriate columns of a sparse matrix
            r = collect(1:m)''  # 91
            r = r[:, fill(1, k-order[ii]+1)]  # 91
            c = collect((order[ii] - k:0) - (order[ii] - minorder))'  # 92
            c = c[fill(1, m), :] + ind[:, fill(1, k-order[ii] + 1)]  # 93
            B[ii] = sparse(vec(r), vec(c), vec(bas[:, 1:k-order[ii]+1]),
                           m, n-order[ii])  # 94 TODO: using vec might be slow

            # 96-100
            if order[ii] > 0
                B[ii] = B[ii] * D[order[ii]]
            elseif order[ii] < 0
                B[ii] = B[ii] * I[-order[ii]]
            end
        end
    end

    # 105
    if length(order) == 1
        B = B[1]
    end
    B, x
end

# ------------------- #
# Lin-family routines #
# ------------------- #

# lindef.m -- DONE
function lindef(breaks::Vector, evennum::Int=0)
    n = length(breaks)  # 28
    !(issorted(breaks)) && error("breaks should be increasing")

    if evennum != 0
        if length(breaks) == 2
            breaks = linspace(breaks[1], breaks[2], evennum)
        else
            if length(breaks) < 2
                error("breaks must have at least 2 elements")
            end

            if any(abs(diff(diff(breaks))) .> 5e-15*mean(abs(breaks)))
                error("Breaks not evenly spaced")
            end
            evennum = length(breaks)
        end
    end
    n = length(breaks)
    a = breaks[1]
    b = breaks[end]
    n, a, b, Any[breaks, evennum]
end

# linnode.m -- DONE
linnode(breaks, evennum) = breaks

# lindop.m
function lindop(breaks, evennum, order)
    # TODO: fill in lindop when I need this
    nothing
end


# linbase.m -- DONE
function linbase(breaks, evennum=0, x=breaks, order=0)
    n = length(breaks)

    # If multiple orders are requested make recursive call 35-44
    k = length(order)
    if k > 1
        B = cell(k)
        for ii=1:k
            B[ii] = linbase(breaks, evennum, x, order[ii])[1]
        end
        return B, x
    end

    # 46-49
    if order != 0
        D, n, a, b, parms = lindop(breaks, evennum, order)
        B = linbase(parms..., x)[1] * D[end]
        return B, x
    end

    m = size(x, 1)

    # Determine the maximum index of
    #   the breakpoints that are less than or equal to x,
    #   (if x=b use the index of the next to last breakpoint).
    if evennum != 0
        ind = fix((x-breaks[1]).*((n-1)./(breaks[end]-breaks[1]))) + 1
        ind = min(max(ind, 1), n-1)
    else
        ind = lookup(breaks, x, 3)
    end

    z = (x-breaks[ind])./(breaks[ind+1]-breaks[ind])
    B = sparse(vcat(1:m, 1:m), vcat(ind, ind+1), vcat(1-z, z), m, n)
    return B, x
end
