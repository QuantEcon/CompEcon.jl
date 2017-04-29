# ------------ #
# Original API #
# ------------ #

# lindef.m -- DONE
function lindef(breaks::AbstractVector, evennum::Int=0)
    p = LinParams(breaks, evennum)
    return (
        length(p.breaks),
        minimum(p.breaks),
        maximum(p.breaks),
        Any[p.breaks, evennum],
        p
    )
end

# linnode.m -- DONE
linnode(breaks, evennum) = nodes(LinParams(breaks, evennum))

# lindop.m
function lindop(breaks, evennum, order)
    D, params = derivative_op(LinParams(breaks, evennum), order)
    n, a, b = length(params.breaks), params.breaks[1], params.breaks[end]
    D, n, a, b, params
end

# linbase.m -- DONE
function linbase(breaks, evennum=0, x=breaks, order=0)
    B = evalbase(LinParams(breaks, evennum), x, order)
    B, x
end
