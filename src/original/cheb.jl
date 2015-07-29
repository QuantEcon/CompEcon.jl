# ------------ #
# Original API #
# ------------ #

# chebdef.m -- DONE
function chebdef(n::Int, a::Real, b::Real)
    cb = Basis(Cheb(), n, a, b)
    return cb.n[1], cb.a[1], cb.b[1], Any[n, a, b]
end


# chebnode.m -- DONE
chebnode(n::Int, a::Real, b::Real, nodetype=0) =
    nodes(ChebParams(n, a, b), nodetype)

# chebdop.m -- DONE
chebdop(n, a, b, order=1) =
    derivative_op(ChebParams(n, a, b), order)


# chebbase.m -- DONE
chebbase(n, a, b, x=chebnode(n, a, b, 1), order=0, nodetype=1) =
    evalbase(ChebParams(n, a, b), x, order, nodetype)

# chebbasex.m -- DONE
chebbasex(n, a, b, x) = evalbasex(ChebParams(n, a, b), x)
