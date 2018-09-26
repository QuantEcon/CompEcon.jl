using CompEcon
using Dierckx


# Create functions to test against
f1(x, y) = sin(x).*cos(y)
f2(x, y) = exp(x) .* cos(y)
f3(x, y) = 1.5 + log(x .* y)
fs = [f1, f2, f3]
nf = length(fs)

x = collect(range(.5, stop=2, length=35))
y = collect(range(.2, stop=4, length=35))

function compare_D_CE_levels(f::Function, x::Array{Float64, 1}, y::Array{Float64, 1}, order::Int)

    # Create Dierckx Spline (and data)
    n = length(x)
    z = f(x, y')
    dspl = Spline2D(x, y, z; s=0., kx=order, ky=order)

    # Create CompEcon Spline (and data)
    cebasis_x = Basis(SplineParams(x, 0, order))
    cebasis_y = Basis(SplineParams(y, 0, order))
    cebasis = Basis(cebasis_x, cebasis_y)
    xx = nodes(cebasis)[1]
    yy = f(xx[:, 1], xx[:, 2])
    c, bs = funfitxy(cebasis, xx, yy)

    # Evaluate Splines on finer grid
    xfine = collect(range(x[1], stop=x[end], length=2*n + 1))
    yfine = collect(range(y[1], stop=y[end], length=2*n + 1))
    zfine = f(xfine, yfine')

    deval = Dierckx.evalgrid(dspl, xfine, yfine)
    derr = abs(zfine - deval)

    ceeval = CompEcon.funeval(c, cebasis, gridmake(xfine, yfine))
    ceerr = abs(zfine - reshape(ceeval, 2*n+1, 2*n+1))

    return derr, ceerr
end


for i=1:nf

    curr_f = fs[i]

    derr1_levels, ceerr1_levels = compare_D_CE_levels(curr_f, x, y, 1)
    derr2_levels, ceerr2_levels = compare_D_CE_levels(curr_f, x, y, 2)
    derr3_levels, ceerr3_levels = compare_D_CE_levels(curr_f, x, y, 3)

    println("Linear Approximation to Function")
    println("Dierckx Level Errors \n\t Max: $(maximum(derr1_levels)) \n\t Min: $(minimum(derr1_levels)) \n\t Mean: $(mean(derr1_levels))")
    println("CompEcon Level Errors \n\t Max: $(maximum(ceerr1_levels)) \n\t Min: $(minimum(ceerr1_levels)) \n\t Mean: $(mean(ceerr1_levels))")

    println("Quadratic Approximation to Function")
    println("Dierckx Level Errors \n\t Max: $(maximum(derr2_levels)) \n\t Min: $(minimum(derr2_levels)) \n\t Mean: $(mean(derr2_levels))")
    println("CompEcon Level Errors \n\t Max: $(maximum(ceerr2_levels)) \n\t Min: $(minimum(ceerr2_levels)) \n\t Mean: $(mean(ceerr2_levels))")

    println("Cubic Approximation to Function")
    println("Dierckx Level Errors \n\t Max: $(maximum(derr3_levels)) \n\t Min: $(minimum(derr3_levels)) \n\t Mean: $(mean(derr3_levels))")
    println("CompEcon Level Errors \n\t Max: $(maximum(ceerr3_levels)) \n\t Min: $(minimum(ceerr3_levels)) \n\t Mean: $(mean(ceerr3_levels))")

end