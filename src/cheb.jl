# ---------- #
# Julian API #
# ---------- #

function Basis(::Cheb, n::Int, a::Real, b::Real)  # chebdef
    n <= 0 && error("n must be positive")
    a >= b && error("left endpoint (a) must be less than right end point (b)")
    a, b = map(Float64, (a, b))
    Basis(Cheb(), n, a, b, ChebParms(n, a, b))
end

# define methods for ChepParms type
Basis(p::ChebParms) = Basis(Cheb(), p.n, p.a, p.b)

nodes(p::ChebParms, nodetype=0) = chebnode(p.n, p.a, p.b, nodetype)

derivative_op(p::ChebParms, order=1) = chebdop(p.n, p.a, p.n, order)

evalbase(p::ChebParms, x=nodes(p, 1), order=0, nodetype=1) =
    chebbase(p.n, p.a, p.b, x, order, nodetype)

evalbasex(p::ChebParms, x) = chebbasex(p.n, p.a, p.b, x)

