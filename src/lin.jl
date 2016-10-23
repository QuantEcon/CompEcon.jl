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
Basis(p::LinParams) = Basis(Lin(), p.breaks, p.evennum)

nodes(p::LinParams) = p.breaks

function derivative_op(p::LinParams, order::Int=1)
    breaks, evennum = p.breaks, p.evennum

    newbreaks = breaks
    n = length(breaks)
    D = Array(SparseMatrixCSC{Float64, Int}, abs(order))

    for i in 1:order
        d = 1./diff(newbreaks)
        d = sparse([1:n-1; 1:n-1], [1:n-1; 2:n], [-d; d], n-1, n)
        if i > 1
            D[i] = d*D[i-1]
        else
            D[1] = d
        end
        newbreaks = (newbreaks[1:end-1]+newbreaks[2:end])/2
        n = n-1
    end

    for i in -1:-1:order
        newbreaks=[dot([3, -1], newbreaks[1:2]);
                   (newbreaks[1:end-1]+newbreaks[2:end]);
                   dot([-1, 3], newbreaks[end-1:end])]/2
        d = diff(newbreaks)'
        n = n+1
        d = tril(repmat(d, n, 1), -1)
        if i<-1
            D[-i] = d*D[-i-1]
        else
            D[1] = d
        end
        #adjustment to make value at original left endpoint equal 0
        if evennum > 0
            temp = evalbase(LinParams(newbreaks, length(newbreaks)),
                            breaks[1], 0)*D[-i]
        else
            temp = evalbase(LinParams(newbreaks, 0), breaks[1], 0)*D[-i]
        end
        D[-i] = D[-i]-repmat(temp, length(newbreaks), 1)
    end

    params = LinParams(newbreaks, evennum > 0 ? evennum : 0)

    return D, params

end

function _prep_evalbase(p::LinParams, x::Real)
    n = length(p.breaks)
    if p.evennum != 0
        _ind = fix((x-p.breaks[1])*((n-1)/(p.breaks[end]-p.breaks[1])))
        ind = clamp(_ind+1, 1, n-1)
    else
        ind = lookup(p.breaks, x, 3)
    end
    1, n, [ind]
end

function _prep_evalbase(p::LinParams, x::AbstractArray)
    m = size(x, 1)
    n = length(p.breaks)

    # Determine the maximum index of
    #   the breakpoints that are less than or equal to x,
    #   (if x=b use the index of the next to last breakpoint).
    if p.evennum != 0
        ind = fix((x-p.breaks[1]).*((n-1)./(p.breaks[end]-p.breaks[1]))) + 1
        clamp!(ind, 1, n-1)
    else
        ind = lookup(p.breaks, x, 3)
    end

    return m, n, ind
end

function evalbase(::Type{SparseMatrixCSC}, p::LinParams,
                  x::Union{Real,AbstractArray}=nodes(p), order::Int=0)
    # 46-49
    if order != 0
        D, params = derivative_op(p, order)
        B = evalbase(params, x, 0) * D[end]
        return B
    end

    m, n, ind = _prep_evalbase(p, x)
    z = similar(x)
    for i in 1:length(x)
        z[i] = (x[i]-p.breaks[ind[i]])/(p.breaks[ind[i]+1]-p.breaks[ind[i]])
    end

    return sparse(vcat(1:m, 1:m), vcat(ind, ind+1), vcat(1-z, z), m, n)
end

function evalbase(::Type{SplineSparse}, p::LinParams,
                  x::Union{Real,AbstractArray}=nodes(p), order::Int=0)
    # 46-49
    if order != 0
        error("derivatives un-supported right now")
    end

    m, n, ind = _prep_evalbase(p, x)

    z = Array(eltype(x), 2*length(x))
    for i in 1:length(x)
        ix = 2i
        z[ix] = (x[i]-p.breaks[ind[i]])/(p.breaks[ind[i]+1]-p.breaks[ind[i]])
        z[ix-1] = 1 - z[ix]
    end

    return SplineSparse(1, 2, n, z, ind)
end

evalbase(p::LinParams, x::Union{Real,AbstractArray}=nodes(p), order::Int=0) =
    evalbase(SparseMatrixCSC, p, x, order)

function evalbase(p::LinParams, x, order::AbstractArray{Int})
    out = Array(SparseMatrixCSC{Float64,Int}, size(order))

    for I in eachindex(order)
        out[I] = evalbase(p, x, order[I])
    end
    return out
end
