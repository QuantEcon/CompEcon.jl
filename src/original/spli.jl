# ------------ #
# Original API #
# ------------ #

# splidef.m -- DONE
function splidef(breaks, evennum=0, k::Int=3)
    sb = Basis(Spline(), breaks, evennum, k)
    return sb.n[1], sb.a[1], sb.b[1], Any[breaks, evennum, k]
end

# splinode.m  -- DONE
splinode(breaks::Vector, evennum::Int, k::Int=3) =
    nodes(SplineParams(breaks, evennum, k))


function splidop(breaks, evennum=0, k=3, order=1)
    D, p = derivative_op(SplineParams(breaks, evennum, k), order)
    n = length(breaks) + k - 1

    D, n-order, breaks[1], breaks[end], Any[breaks, evennum, k-order]
end

# splibas.m -- DONE
function splibase(breaks::Vector, evennum, k=3, x=splinode(breaks, evennum, k),
                  order=0)
    B = evalbase(SplineParams(breaks, evennum, k), x, order)
    B, x
end
