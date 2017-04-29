# ------------ #
# Original API #
# ------------ #

# splidef.m -- DONE
function splidef(breaks, evennum=0, k::Int=3)
    p = SplineParams(breaks, evennum, k)
    return (
        length(p.breaks),
        minimum(p.breaks),
        maximum(p.breaks),
        Any[p.breaks, 0, k],
        p
    )
end

# splinode.m  -- DONE
splinode(breaks::AbstractVector, evennum::Int, k::Int=3) =
    nodes(SplineParams(breaks, evennum, k))

function splidop(breaks, evennum=0, k=3, order=1)
    D, p = derivative_op(SplineParams(breaks, evennum, k), order)
    n = length(breaks) + k - 1

    D, n-order, breaks[1], breaks[end], Any[breaks, evennum, k-order]
end

# splibas.m -- DONE
function splibase(breaks::AbstractVector, evennum, k=3, x=splinode(breaks, evennum, k),
                  order=0)
    B = evalbase(SplineParams(breaks, evennum, k), x, order)
    B, x
end
