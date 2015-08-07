# --------------- #
# Chebyshev Basis #
# --------------- #

function Basis(::Cheb, n::Int, a::Real, b::Real)  # chebdef
    n <= 0 && error("n must be positive")
    a >= b && error("left endpoint (a) must be less than right end point (b)")
    Basis(Cheb(), n, a, b, ChebParams(n, a, b))
end

# define methods for ChepParams type
Basis(p::ChebParams) = Basis(Cheb(), p.n, p.a, p.b)

# chebnode.m -- DONE
function nodes(p::ChebParams, nodetype=0)
    n, a, b = p.n, p.a, p.b
    s = (b-a) / 2  # 21
    m = (b+a) / 2  # 22

    if nodetype < 2  # usual nodes
        k = π*(0.5:(maximum(n) - 0.5))'  # 25
        x = m - cos(k[1:n] / n) * s  # 26
        if nodetype == 1  # Extend nodes to endpoints
            aa = x[1]  # 28
            bb = x[end]  # 29
            x = (bb*a - aa*b)/(bb-aa) + (b-a)/(bb-aa)*x  # 30
        end
    else # Lobatto nodes
        k = pi*(0:(maximum(n)-1))  # 33
        x = m - cos(k[1:n] / (n-1)) * s  # 34
    end
    x
end

function derivative_op(p::ChebParams, order=1)
    n, a, b = p.n, p.a, p.b
    if order > 0
        # TODO: figure out some caching mechanism that avoids globals
        D = cell(max(2, order), 1)  # 49
        i = repeat(collect(1:n), outer=[1, n]); j = i'  # 50

        # 51
        inds = find((rem(i + j, 2) .== 1) & (j.>i))
        r, c = similar(inds), similar(inds)
        for i in 1:length(inds)
            r[i], c[i] = ind2sub((n, n), inds[i])
        end

        d = sparse(r, c, (4/(b-a)) * (vec(j[1, c])-1), n-1, n)  # 52
        d[1, :] ./= 2  # 53
        D[1] = d  # 54
        for ii=2:max(2, order)
            D[ii] = d[1:n-ii, 1:n-ii+1] * D[ii-1]  # 56
        end
    else
        D = cell(abs(order), 1)  # 64
        nn = n - order  # 65
        i = (0.25 * (b - a)) ./(1:nn)  # 66
        d = sparse(vcat(1:nn, 1:nn-2), vcat(1:nn, 3:nn), vcat(i, -i[1:nn-2]),
                   nn, nn)  # 67
        d[1, 1] *= 2  # 68
        d0 = ((-1).^(0:nn-1)') .* sum(d, 1)  # 69
        D[1] = sparse(vcat(d0[1:n]', d[1:n, 1:n]))  # 70
        for ii=-2:-1:order
            ind = 1:n-ii-1
            D[-ii] = sparse(vcat(d0[ind]', d[ind, ind]) * D[-ii-1])
        end
    end
    D, n-order, a, b, Any[n, a, b]
end

function evalbase(p::ChebParams, x=nodes(p, 1), order=0, nodetype=1)
    n, a, b = p.n, p.a, p.b
    minorder = min(0, minimum(order))  # 30

    # compute 0-order basis
    if nodetype == 0
        temp = ((n-0.5):-1:0.5)''  # 41
        bas = cos((pi/n)*temp.*(0:(n-1-minorder))')  # 42
    else
        bas = chebbasex(n-minorder, a, b, x)  # 44
    end

    if length(order) == 1
        if order != 0
            D = chebdop(n, a, b, order)[1]
            B = bas[:, 1:n-order]*D[abs(order)]
        else
            B = bas
        end
    else
        B = cell(length(order))
        maxorder = maximum(order)
        if maxorder > 0 D = chebdop(n, a, b, maxorder)[1] end
        if maxorder < 0 I = chebdop(n, a, b, minorder)[1] end
        for ii=1:length(order)
            if order[ii] == 0
                B[ii] = bas[:, 1:n]
            elseif order[ii] > 0
                B[ii] = bas[:, 1:n-order[ii]] * D[order[ii]]
            else
                B[ii] = bas[:, 1:n-order[ii]] * I[-order[ii]]
            end
        end
    end

    return B, x
end

function evalbasex(p::ChebParams, x)
    n, a, b = p.n, p.a, p.b
    z = (2/(b-a)) * (x-(a+b)/2)
    m = size(z, 1)
    bas = Array(Float64, m, n);
    bas[:, 1] = 1.0
    bas[:, 2] = z
    z = 2 * z
    @inbounds for j=3:n
        for i=1:m
            bas[i, j] = z[i] .* bas[i, j-1] - bas[i, j-2]
        end
    end
    bas
end
