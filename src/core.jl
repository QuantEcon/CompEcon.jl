# -------------------- #
# Dispatcher functions #
# -------------------- #

const BASE_TYPES = [:spli, :cheb, :lin]
const ABSR_MAP = Dict(
    :none => Direct(),
    :direct => Direct(),
    :tensor => Tensor(),
    :expanded => Expanded(),
)
get_bformat(b::T) where T<:BasisMatrix{Direct} = :direct
get_bformat(b::T) where T<:BasisMatrix{Expanded} = :expanded
get_bformat(b::T) where T<:BasisMatrix{Tensor} = :tensor

function to_dict(bm::BasisMatrix)
    B = Dict{Symbol, Any}()
    B[:order] = bm.order
    B[:format] = get_bformat(bm)
    B[:vals] = bm.vals
    B
end

function bm_from_dict(B::Dict)
    arr_type = eltype(B[:vals])
    bm = BasisMatrix{typeof(ABSR_MAP[B[:format]]),arr_type}(B[:order], B[:vals])
    bm
end

base_exists(s::Symbol) = s in BASE_TYPES

basedef(s::Symbol, args...) =
    s == :spli ? splidef(args...) :
    s == :cheb ? chebdef(args...) :
    s == :lin  ? lindef(args...)  :
    error("somehow you snuck through here you 👺")

basenode(s::Symbol, args...) =
    s == :spli ? splinode(args...) :
    s == :cheb ? chebnode(args...) :
    s == :lin  ? linnode(args...)  :
    error("somehow you snuck through here you 👺")

BasisMatrices.evalbase(s::Symbol, args...) =
    s == :spli ? splibase(args...) :
    s == :cheb ? chebbase(args...) :
    s == :lin  ? linbase(args...)  :
    error("somehow you snuck through here you 👺")

# Helper function
function squeeze_trail(x::AbstractArray)
    sz = size(x)
    squeezers = Int[]
    n = length(sz)
    for i=n:-1:1
        if sz[i] == 1
            push!(squeezers, i)
        else
            break
        end
    end
    squeeze(x, tuple(squeezers...))
end


# ---------------------------- #
# Generic translated functions #
# ---------------------------- #

# from fundef.m  -- DONE
function fundef(foo...)
    d = length(foo)  # 89
    n = zeros(Int, d)  # 93
    b = zeros(d)  # 94
    a = zeros(d)  # 95
    p = Array{Any}(undef, d)  # 96
    _params = Array{BasisMatrices.BasisParams}(undef, d)

    basetype = Array{Symbol}(undef, d)
    for j=1:d
        basetype[j] = foo[j][1]  # 99
        !(base_exists(basetype[j])) && error("Unknown basis $(foo[j][1])")
        n[j], a[j], b[j], p[j], _params[j] = basedef(basetype[j], foo[j][2:end]...)  # 124
    end

    # package output. Lines 143-150
    g = Dict{Symbol, Any}()
    g[:d] = d
    g[:n] = n
    g[:a] = a
    g[:b] = b
    g[:basetype] = basetype
    g[:params] = p
    g[:_basis_params] = _params
    g[:_basis] = Basis(_params...)
    g
end

# fundefn.m
function fundefn(basistype::Symbol, n, a, b, order=3)
    d = length(n)
    length(a) != d && error("a must be same dimension as n")
    length(b) != d && error("b must be same dimension as n")
    any(a .> b) && error("left endpoints must be less than right endpoints")
    any(n .< 2) && error("n(i) must be greater than 1")

    params = Array{Any}(undef, 1, d)
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
funnode(basis::Dict) = nodes(basis[:_basis])

# funbase.m -- DONE
function funbase(basis::Dict, x=funnode(basis)[1], order=fill(0, 1, basis[:d]))
    BasisMatrix(basis[:_basis], Expanded(), x, order).vals[1]
end

# funbasex.m -- DONE
function funbasex(basis::Dict{Symbol}, x=funnode(basis)[1], order=0,
                  bformat::Symbol=:none)
    to_dict(BasisMatrix(basis[:_basis], ABSR_MAP[bformat], x, order))
end

# funfitf.m -- DONE
funfitf(basis, f::Function, args...) = funfitf(basis[:_basis], f, args...)

# funfitxy.m -- DONE
function funfitxy(basis, x, y)
    c, bm = funfitxy(basis[:_basis], x, y)
    return c, to_dict(bm)
end

# funeval.m -- DONE
function funeval(c, basis::Dict, B, _order=0)
    isempty(c) && error("missing basis coefficients")
    order = BasisMatrices._check_order(basis[:d], _order)

    if isa(B, Dict)  # B is a basis structure
        bm = bm_from_dict(B)
        y = funeval(c, bm, order)
        return y, B
    else
        bm = BasisMatrix(basis[:_basis], B, order)
        y = funeval(c, bm, bm.order)
        return y, to_dict(bm)
    end
end

# fund.m
function fund(c, basis, x, hess_opt)
    # TODO: come back when I need this. I think I should probably do something
    #       like what optim does and write functions `f`, `fg!` and `fgh!` to
    #       replicate this for the type of basis instead of this function
    nothing
end

# funbconv.m  -- DONE
function funbconv(b::Dict, order=fill(0, 1, size(b[:order], 2)),
                  format::Symbol=:expanded)
    bm = bm_from_dict(b)
    new_bm = convert(typeof(ABSR_MAP[format]), bm, order)
    to_dict(new_bm)
end
