# ------------ #
# Original API #
# ------------ #

# lindef.m -- DONE
function lindef(breaks::Vector, evennum::Int=0)
    lb = Basis(Lin(), breaks, evennum)
    return lb.n[1], lb.a[1], lb.b[1], Any[breaks, evennum]
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
        D, n, a, b, params = lindop(breaks, evennum, order)
        B = linbase(params..., x)[1] * D[end]
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
