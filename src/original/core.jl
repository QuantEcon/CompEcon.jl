# -------------------- #
# Dispatcher functions #
# -------------------- #

const BASE_TYPES = [:spli, :cheb, :lin]

base_exists(s::Symbol) = s in BASE_TYPES


basedef(s::Symbol, args...) =
    s == :spli ? splidef(args...) :
    s == :cheb ? chebdef(args...) :
    s == :lin  ? lindef(args...)  :
    error("somehow you snuck through here you :japanese_goblin:")

basenode(s::Symbol, args...) =
    s == :spli ? splinode(args...) :
    s == :cheb ? chebnode(args...) :
    s == :lin  ? linnode(args...)  :
    error("somehow you snuck through here you :japanese_goblin:")

evalbase(s::Symbol, args...) =
    s == :spli ? splibase(args...) :
    s == :cheb ? chebbase(args...) :
    s == :lin  ? linbase(args...)  :
    error("somehow you snuck through here you :japanese_goblin:")


# ---------------------------- #
# Generic translated functions #
# ---------------------------- #

# from fundef.m  -- DONE
function fundef(foo...)
    d = length(foo)  # 89
    n = zeros(Int, d)  # 93
    b = zeros(d)  # 94
    a = zeros(d)  # 95
    p = cell(d)  # 96

    basetype = Array(Symbol, d)
    for j=1:d
        basetype[j] = foo[j][1]  # 99
        !(base_exists(basetype[j])) && error("Unknown basis $(foo[j][1])")
        n[j], a[j], b[j], p[j] = basedef(basetype[j], foo[j][2:end]...)  # 124
    end

    # package output. Lines 143-150
    g = Dict{Symbol, Any}()
    g[:d] = d
    g[:n] = n
    g[:a] = a
    g[:b] = b
    g[:basetype] = basetype
    g[:parms] = p
    g
end


# fundefn.m
function fundefn(basistype::Symbol, n, a, b, order=3)
    d = length(n)
    length(a) != d && error("a must be same dimension as n")
    length(b) != d && error("b must be same dimension as n")
    any(a .> b) && error("left endpoints must be less than right endpoints")
    any(n .< 2) && error("n(i) must be greater than 1")

    params = cell(1, d)
    if basistype == :cheb
        for i=1:d params[i] = Any[:cheb, n[i], a[i], b[i]] end
    elseif basistype == :spli
        for i=1:d params[i] = Any[:spli, [a[i], b[i]], n[i]-order+1, order] end
    elseif basistype == :lin
        for i=1:d params[i] = Any[:lin, [a[i], b[i]], n[i]] end
    end

    fundef(params...)
end

# funnode.m -- DONE
function funnode(basis::Dict)
    d = basis[:d]  # 18
    if d == 1  # 19
        xcoord = basenode(basis[:basetype][1], basis[:parms][1]...)  # 20
        x = xcoord  # 21
    else  # 22
        xcoord = cell(d)  # 23
        for j=1:d  # 24
            xcoord[j] = basenode(basis[:basetype][j], basis[:parms][j]...)  #25
        end  # 26
        x = gridmake(xcoord...)  # 27
    end

    return x, xcoord
end

# funbase.m -- DONE
funbase(basis::Dict, x=funnode(basis)[1],
        order=fill(0, 1, basis[:d])) =
    funbasex(basis, x, order, :expanded)[:vals][1]


# funbasex.m -- DONE
function funbasex(basis, x=funnode(basis)[1], order=0, bformat::Symbol=:none)
    d = length(basis[:n])  # 59
    if d > 1 && size(order, 2) == 1  # 62
        order = order * ones(Int, 1, d)  # 63
    end  # 64

    # initialize basis structure (66-74)
    m = size(order, 1)
    if m > 1
        minorder = fill(minimum(order), 1, d)
        numbases = (maximum(order) - minorder) + 1
    else
        minorder = order + zeros(Int, 1, d)
        numbases = fill(1, 1, d)
    end

    B = Dict{Symbol, Any}()  # 76
    B[:order] = minorder
    B[:format] = bformat
    B[:vals] = cell(maximum(numbases), d)  # 77

    # 83-89
    if bformat == :none
        if isa(x, Array{Any})
            bformat = :tensor
        else
            bformat = :direct
        end
    end

    # 91-101
    if d > 1
        if !(isa(x, Array{Any})) && bformat == :tensor
            error("Must pass a cell array to form a tensor format basis")
        end
        if isa(x, Array{Any}) && bformat == :direct
            # it would be more efficient in this case to
            # use the tensor form to compute the bases and then
            # to use indexing to expand to the direct form
            #      x=gridmake(x); % convert to grid for direct form
        end
    end

    #103-107
    B[:format] = is(x, Array{Any}) ? :tensor : :direct

    if B[:format] == :tensor  # 111
        for j=1:d
            # 113-117
            if (m > 1)
                orderj = unique(order[:, j])
            else
                orderj = order[1, j]
            end

            #118-122
            if length(orderj) == 1
                B[:vals][1, j] = evalbase(basis[:basetype][j],
                                          basis[:parms][j]..., x[j], orderj)
            else
                B[:vals][orderj-minorder[j]+1, j] =
                    evalbase(basis[:basetype][j], basis[:parms][j]..., x[j],
                             orderj)
            end
        end
    else  # B[:format] == :direct
        for j=1:d
            # 126-130
            if (m > 1)
                orderj = unique(order[:, j])
            else
                orderj = order[1, j]
            end

            #131-135
            if length(orderj) == 1
                B[:vals][1, j] = evalbase(basis[:basetype][j],
                                          basis[:parms][j]..., x[:, j], orderj)[1]
            else
                B[:vals][orderj-minorder[j]+1, j] =
                    evalbase(basis[:basetype][j], basis[:parms][j]..., x[:, j],
                             orderj)[1]
            end
        end
    end

    # if d == 1, switch to expanded format
    # 140-145
    if size(B[:vals], 2) == 1  # 1 dimension
        B[:format] = :expanded
        B[:order] = order
        B[:vals] = B[:vals][collect(order + (1 - minimum(order)))]
        return B
    end

    # create expanded format  148-153
    if bformat == :expanded
        B = funbconv(B, order, :expanded)
    elseif bformat == :direct
        if isa(x, Array{Any})
            B = funbconv(B, order, :direct)
        end
    end
    return B
end

# lookup.m -- DONE
function lookup(tabvals::Vector, x::Vector, endadj=0)
    n = prod(size(x))
    m = length(tabvals)
    if endadj >= 2
        m = m - sum(tabvals .== tabvals[end])
    end

    ind = sortperm(vcat(tabvals[1:m], x))
    temp = find(ind .>m)
    j = ind[temp] - m
    ind = reshape(temp .- (1:n), size(x)...)
    ind[j] = ind[:]

    if endadj == 1 || endadj == 3
        ind[ind .== 0] = sum(tabvals .== tabvals[1])
    end
    ind
end

# funfitf.m -- DONE
function funfitf(basis, f::Function, args...)
    x = funnode(basis)[1]
    y = f(x, args...)
    funfitxy(basis, x, y)[1]
end

# funfitxy.m -- DONE
function funfitxy(basis, x, y)
    m = size(y, 1)
    prod(basis[:n]) > m && error("Cannot be more basis functions than points")

    if isa(x, Dict)  # 61 : use precomputed basis structure
        B = x

        if B[:format] == :tensor
            if any(B[:order][1, :] != 0)
                m = "invalid basis structure - first elements must be order 0"
                error(m)
            end
            B[:vals] = B[:vals][1, :]  # 68
            c = ckronxi(B[:vals], y, basis[:d]:-1:1)  # 66

        elseif B[:format] == :direct
            B = funbconv(B, zeros(1, size(B[:vals], 2)), :expanded)  # 71
            c = B[:vals][1] \ y  # 72

        elseif B[:format] == :expanded
            c = B[:vals][1] \ y  # 74
        end

    elseif isa(x, Array{Any})  # evaluate at grid points
        # 77
        mm = 1
        for i=1:size(x, 2)
            mm = mm * size(x[i], 1)
        end

        mm != m && error("x and y are incompatible")

        B = funbasex(basis, x, 0)  # 81
        c = ckronxi(B[:vals], y, basis[:d]:-1:1)  # 82

    else  # evaluate at arbitrary points
        size(x, 1) != m && error("x and y are incompatible")

        B = funbasex(basis, x, 0, :expanded)
        c = B[:vals][1] \ y
    end
    return c, B
end


# funeval.m -- DONE
function funeval(c, basis, B, order=0)
    isempty(c) && error("missing basis coefficients")

    if !isa(B, Dict)  # B is not a basis structure, so construct one
        d = basis[:d]  # 60
        if size(B, 2) != d  # TODO: this will be ill-defiend for non array B's
            error("x must have d=$(d) columns")  # 62
        end

        if size(order, 2) == 1
            order = order * fill(1, 1, d)  # 65
        end
        B = funbasex(basis, B, order)  # 67
    else
        if size(order, 2) == 1
            order = order * ones(Int, 1, size(B[:order], 2))  # 69
        end
    end

    # now we we have a basis structure
    if isa(B[:vals], Array{Any})
        y = B[:format] == :tensor   ? funeval1(c, B, order) :  # 75
            B[:format] == :direct   ? funeval2(c, B, order) :  # 77
            B[:format] == :expanded ? funeval3(c, B, order) :  # 79
            error("Basis structure has invalid ``format`` field")  # 81
    else
        y = B[:vals] * c  # 84
    end

    return y, B
end

# inside funeval.m. Evaluates according to ``:tensor`` format -- DONE
function funeval1(c, B, order=fill(0, 1, size(B[:order], 2)))
    kk, d = size(order)  # 95
    # 98 reverse the order of evaluation: B(d)xB(d-1)x...xB(1)
    order = flipdim(order+1*(size(B[:vals], 1) * (0:d-1)' - B[:order]+1), 2)

    # 99
    nx = 1
    for j=1:d
        nx = nx * size(B[:vals][i, j], 1)
    end

    f = zeros(nx, size(c, 2), kk)  # 100

    for i=1:kk
        f[:, :, i] = ckronx(B.vals, c, order[i, :])  # 102
    end
    return squeeze_trail(f)
end

# inside funeval.m. Evaluates according to ``:direct`` format -- DONE
function funeval2(c, B, order=fill(0, 1, size(B[:order], 2)))
    kk, d = size(order)  # 95
    # 114 reverse the order of evaluation: B(d)xB(d-1)x...xB(1)
    order = flipdim(order+1*(size(B[:vals], 1) * (0:d-1)' - B[:order]+1), 2)

    f = zeros(size(B[:vals][1], 1), size(c, 2), kk)  # 116

    for i=1:kk
        f[:, :, i] = cdprodx(B[:vals], c, order[i, :])  # 118
    end
    return squeeze_trail(f)
end

# inside funeval.m. Evaluates according to ``:expanded`` format -- DONE
function funeval3(c, B, order=fill(0, 0, 0))
    if isempty(order)
        if isa(B[:vals], Array{Any})
            kk = length(B[:vals])  # 126
            order = 1:kk'  # 127
        else
            kk = 1  # 128
            order = 1  # 129
        end
    else
        kk = size(order, 1)  # 133
    end

    if isa(B[:vals], Array{Any})
        nx = size(B[:vals][1], 1)
        f = zeros(nx, size(c, 2), kk)
        for i=1:kk
            # 140 determine which element of B[:vals] is the desired basis
            ii = Int[]
            for row=1:size(B[:order], 1)
                r = B[:order][row, :]
                if r == order[i, :]
                    push!(ii, row)
                end
            end

            # 141-143
            isempty(ii)  && error("Requested basis matrix is not available")

            length(ii) > 1 &&  warn("redundant request in funeval3")  # 145

            # NOTE: must do even when length[i] == 1 b/c want element of cell
            #       and indexing cell with vector in julia gives cell instead
            #       of the element
            ii = ii[1]  # 146

            f[:, :, i] = B[:vals][ii]*c  #148
         end
     else
        nx = size(B[:vals], 1)  # 151
        f = zeros(nx, size(c, 2), kk)  # 152
        for i=1:kk
            f[:, :, i] = B[:vals]*c  # 154
        end
    end

    return squeeze_trail(f)
end

# fund.m
function fund(c, basis, x, hess_opt)
    # TODO: come back when I need this. I think I should probably do something
    #       like what optim does and write functions `f`, `fg!` and `fgh!` to
    #       replicate this for the type of basis instead of this function
    nothing
end




# funbconv.m  -- DONE
function funbconv(b, order=fill(0, 1, size(b[:order], 2)),
                  format::Symbol=:expanded)

    d = size(b[:order], 2)  # 24
    numbas, d1 = size(order, 1), size(order, 2)  # 34
    d1 != d && error("ORDER incompatible with basis functions")  # 35-37

    # 39-41
    if any(minimum(order) .< b[:order])
        error("Order of derivative operators exceeds basis")
    end

    B = Dict{Symbol, Any}()
    if format == :expanded
        B[:vals] = cell(numbas)
        B[:order] = order
        B[:format] = format

        if b[:format] == :tensor
            m, n = 1, 1
            for j=1:d  # 50
                m = m * size(b[:vals][1, j], 1)  # 51
                n = n * size(b[:vals][1, j], 2)  # 52
            end

            for i=1:numbas  # 54
                B[:vals][i] = b[:vals][order[i, d] - b[:order][d]+1, d]  # 55
                for j=d-1:-1:1  # 56
                    # 57
                    B[:vals][i] = kron(B[:vals][i], b[:vals][order[i, j] - b[:order][j]+1, j])
                end
            end
        elseif b[:format] == :direct
            n = 1  # 61
            for j=1:d  # 61
                n = n * size(b[:vals][1, j], 2)  # 61
            end

            for i=1:numbas
                B[:vals][i] = b[:vals][order[i, d] - b[:order][d]+1, d]  # 63
                for j=d-1:-1:1
                    B[:vals][i] = row_kron(B[:vals][i],
                                           b[:vals][order[i, j] - b[:order][j]+1, j])  #65
                end
            end
        elseif b[:format] == :expanded
            warn("Basis already in expanded form")
            B = b
        else
            error("Improper basis format")
        end

    elseif format == :direct && b[:format] == :tensor
        B[:vals] = cell(numbas, d)  # 75
        B[:order] = order  # 76
        B[:format] = format  # 77
        ind = cell(1, d)  # 78

        for j=1:d, i=1:size(b[:vals], 1)
            if !(isempty(b[:vals][i, j]))  # 81
                ind[j] = collect(1:size(b[:vals][i, j], 1))  # 82
                break
            end
        end

        ind = gridmake(ind...)  # 87

        for j=1:d, i=1:numbas
            if !(isempty(b[:vals][i, j]))  # 90
                B[:vals][i, j] = B[:vals][i, j][ind[j], :]
            end
        end
    else
        error("Not implemented for this option")
    end

    return B
end
