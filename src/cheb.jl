# ------------ #
# Original API #
# ------------ #

# chebdef.m -- DONE
function chebdef(n::Int, a::Real, b::Real)
    p = ChebParams(n, a, b)
    return (p.n, p.a, p.b, Any[p.n, p.a, p.b], p)
end


# chebnode.m -- DONE
chebnode(n::Int, a::Real, b::Real, nodetype=0) =
    nodes(ChebParams(n, a, b), nodetype)

# chebdop.m -- DONE
function chebdop(n, a, b, order=1)
    D, p = derivative_op(ChebParams(n, a, b), order)
    D, p.n-order, p.a, p.b, Any[p.n-order, p.a, p.b]
end


# chebbase.m -- DONE
function chebbase(n, a, b, x=chebnode(n, a, b, 1), order=0, nodetype=1)
    B = evalbase(ChebParams(n, a, b), x, order, nodetype)
    B, x
end

# chebbasex.m -- DONE
chebbasex(n, a, b, x) = evalbasex(ChebParams(n, a, b), x)
