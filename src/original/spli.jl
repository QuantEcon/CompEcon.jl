# ------------ #
# Original API #
# ------------ #

# splidef.m -- DONE
function splidef(breaks, evennum=0, k::Int=3)
    sb = Basis(Spline(), breaks, evennum, k)
    return sb.n[1], sb.a[1], sb.b[1], Any[breaks, evennum, k]
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
