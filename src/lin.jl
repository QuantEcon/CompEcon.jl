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

nodes(p::LinParams) = linnode(p.breaks, p.evennum)

derivative_op(p::LinParams, order=1) = lindop(p.breaks, p.evennum, order)

evalbase(p::LinParams, x=nodes(p), order=0) =
    linbase(p.breaks, p.evennum, x, order)
