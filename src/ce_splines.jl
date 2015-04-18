# TODO: no cheb functions have been written yet.

const BASE_TYPES = [:spli, :cheb]

base_exists(s::Symbol) = s in BASE_TYPES

function basedef(s::Symbol, args...)
    if s == :spli
        return splidef(args...)
    elseif s == :cheb
        return chebdef(args...)
    else
        error("somehow you snuck through here you :japanese_goblin:")
    end
end


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
            D = chebdop(n, a, b, order)
            B = bas[:, 1:n-order]*D[abs(order)]
        else
            B = bas
        end
    else
        B = cell(length(order))
        maxorder = maximum(order)
        if maxorder > 0 D = chebdop(n, a, b, maxorder) end
        if maxorder < 0 I = chebdop(n, a, b, minorder) end
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

    return B
end

# chebdop.m
function chebdop(n, a, b, order=1)
    if order > 0
        # Use previously stored values for order=1 and order=2
        # this speeds up calculations considerably with repeated calls

        # TODO: come back here when I need this.
    end
    D, n, a, b, parms
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


# from fundef.m  -- DONE
function fundef(foo...)
    d = length(foo)  # 89
    n = zeros(d)  # 93
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

function basenode(s::Symbol, args...)
    if s == :spli
        return splinode(args...)
    elseif s == :cheb
        return chebnode(args...)
    else
        error("somehow you snuck through here you :japanese_goblin:")
    end
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
        x = gridmake(xcoord)  # 27
    end

    return x, xcoord
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

function evalbase(s::Symbol, args...)
    if s == :spli
        return splibase(args...)
    elseif s == :cheb
        return chebbase(args...)
    else
        error("somehow you snuck through here you :japanese_goblin:")
    end
end


# funbasex.m -- DONE
function funbasex(basis, x=funnode(basis)[1], order=0, bformat::Symbol=:none)
    d = length(basis[:n])  # 59
    if d > 1 && size(order, 2) == 1  # 62
        order = order * ones(1, d)  # 63
    end  # 64

    # initialize basis structure (66-74)
    m = size(order, 1)
    if m > 1
        minorder = fill(minimum(order), d)
        numbases = (maximum(order) - minorder) + 1
    else
        minorder = order + zeros(d)
        numbases = fill(1, d)
    end

    B = Dict{Symbol, Any}()  # 76
    B[:vals] = cell(maximum(numbases), d)  # 77

    # 83-89
    if bformat == :none
        if isa(x, Vector{Any})
            bformat = :tensor
        else
            bformat = :direct
        end
    end

    # 91-101
    if d > 1
        if !(isa(x, Vector{Any})) && bformat == :tensor
            error("Must pass a cell array to form a tensor format basis")
        end
        if isa(x, Vector{Any}) && bformat == :direct
            # it would be more efficient in this case to
            # use the tensor form to compute the bases and then
            # to use indexing to expand to the direct form
            #      x=gridmake(x); % convert to grid for direct form
        end
    end

    #103-107
    if is(x, Vector{Any})
        B[:format] = :tensor
    else
        B[:format] = :direct
    end

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
        B[:vals] = B[:vals][order + (1 - minimum(order))]
        return B
    end

    # create expanded format  148-153
    if bformat == :expanded
        B = funbconv(B, order, :expanded)
    elseif bformat == :direct
        if isa(x, Vector{Any})
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





#=
function funbconv(b, order=fill(0, size(b[:order], 2)),
                  format::Symbol=:expanded)

    d = size(b[:order], 2)  # 24
    numbas, d1 = size(order, 1), size(order, 2)  # 34
    d1 != d && error("ORDER incompatible with basis functions")  # 35-37

    # 39-41
    if any(minimum(order) .< b.order)
        error("Order of derivative operators exceeds basis")
    end

    B = Dict{Symbol, Any}()
    if format == :expanded
        B[:vals] = cell(numbas)
        B[:order] = order
        B[:format] = format

        if b[:format] == :tensor
            n, n = 1, 1
            for j=1:d
                m = m * size(b[:vals][1])
                # TODO: pick up here when I need this

end
=#


# splidop.m
function splidop(breaks, evennum=0, k=3, order=1)
    any(order .>= k) && error("Order of differentiation must be less than k")

    # 38-40
    n = length(breaks) + k - 1
    kk = max(k - 1, k - order - 1)
    augbreaks = vcat(fill(breaks[1], kk), breaks, fill(breaks[end], kk))

    D = cell(abs(order))
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



# funbase.m
function funbase(basis::Dict, x=funnode(basis)[1],
                 order::Vector{Int}=zeros(Int, basis[:d]))
    B = funbasex(basis, x, order, :expanded)
    B[:vals][1]
end


# TODO: still need to write funeval, funfitf, finfitxy, fund
# TODO: also need `lin` family: linbase, lindef, lindop, linnode, minterp
