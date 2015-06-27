# ---------------------- #
# Piecewise linear Basis #
# ---------------------- #

function Basis(::Lin, breaks::Vector, evennum::Int=0)
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
    Basis(Lin(), n, a, b, LinParams(breaks, evennum))
end

# define methods for LinParams type
Basis(p::LinParams) = Basis(Spline(), p.breaks, p.evennum)

nodes(p::LinParams) = p.breaks

function derivative_op(p::LinParams, order=1)
    # TODO: fill in lindop when I need this
    nothing
end

function evalbase(p::LinParams, x=nodes(p), order=0)
    breaks, evennum = p.breaks, p.evennum
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
