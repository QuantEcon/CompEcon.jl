# ------------ #
# Original API #
# ------------ #

# lindef.m -- DONE
function lindef(breaks::Vector, evennum::Int=0)
    lb = Basis(Lin(), breaks, evennum)
    return lb.n[1], lb.a[1], lb.b[1], Any[breaks, evennum]
end

# linnode.m -- DONE
linnode(breaks, evennum) = nodes(LinParam(breaks, evennum))

# lindop.m
lindop(breaks, evennum, order) =
    derivative_op(LinParam(breaks, evennum), order)

# linbase.m -- DONE
linbase(breaks, evennum=0, x=breaks, order=0) =
    evalbase(LinParam(breaks, evennum), x, order)
