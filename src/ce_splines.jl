# TODO: still need to write fund, minterp
# TODO: also need splidop, lindop
# TODO: funeval fails for scalar input and does weird thing for 1-element
#       vector input

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
