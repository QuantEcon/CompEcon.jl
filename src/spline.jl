# ---------- #
# Julian API #
# ---------- #

include("original/spli.jl")

function Basis(::Spline, breaks::Vector, evennum::Int, k::Int)
    # error handling
    k < 0 && error("spline order must be positive")
    length(breaks) < 2 && error("Must have at least two breakpoints")
    any(diff(breaks) .< 0) && error("Breakpoints must be non-decreasing")

    if evennum == 0  # 43
        if length(breaks) == 2  # 44
            evennum = 2
        end
    else
        if length(breaks) == 2
            breaks = linspace(breaks[1], breaks[2], evennum)
        else
            error("Breakpoint squence must contain 2 values when evennum > 0")
        end
    end

    n = length(breaks) + k - 1
    a = breaks[1]
    b = breaks[end]
    Basis(Spline(), n, a, b, SplineParms(breaks, evennum, k))
end

Basis(p::SplineParms) = Basis(Spline(), p.breaks, p.evennum, p.k)

nodes(p::SplineParms) = splinode(p.breaks, p.evennum, p.k)
