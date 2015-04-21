function golden(f::Function, a::AbstractVector, b::AbstractVector;
                tol=eps()/10, maxit=1000)

    α1 = (3 - sqrt(5)) / 2
    α2 = 1 - α1
    d = b - a
    x1 = a + α1*d
    x2 = a + α2*d
    s = ones(a)
    f1 = f(x1)
    f2 = f(x2)

    d = α1*α2*d
    it = 0

    while any(d .> tol) && it < maxit
        it += 1
        i = f2 .> f1
        x1[i] = x2[i]
        f1[i] = f2[i]
        d *= α2
        x2 = x1 + s.*(i- ~i).*d
        s = sign(x2 - x1)
        f2 = f(x2)
    end

    it >= maxit && warn("`golden`: maximum iterations exceeded")

    i = f2 .> f1
    x1[i] = x2[i]
    f1[i] = f2[i]

    x1, f1
end

function golden(f::Function, a::Real, b::Real; tol=eps()/10, maxit=1000)
    x1, f1 = golden(f, [a], [b]; tol=tol, maxit=maxit)
    x1[1], f1[1]
end
