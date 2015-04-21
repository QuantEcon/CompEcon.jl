# ---------- #
# Julian API #
# ---------- #

include("original/cheb.jl")

function Basis(::Cheb, n::Int, a::Real, b::Real)  # chebdef
    n <= 0 && error("n must be positive")
    a >= b && error("left endpoint (a) must be less than right end point (b)")
    a, b = map(Float64, (a, b))
    Basis(Cheb(), n, a, b, ChebParms(n, a, b))
end

Basis(p::ChebParms) = Basis(Cheb(), p.n, p.a, p.b)

nodes(p::ChebParms, nodetype=0) = chebnode(p.n, p.a, p.b, nodetype)
