using CompEcon
using Dierckx


# Create functions to test against
f1(x) = sin(x)
fp1(x) = cos(x)
f2(x) = exp(x)
fp2(x) = exp(x)
f3(x) = 1.5 + log(x)
fp3(x) = 1./x
fs = [f1, f2, f3]
fps = [fp1, fp2, fp3]
nf = length(fs)

x = collect(linspace(.5, 2, 35))
finex = collect(linspace(1e-2, 2*pi - 1e-2, 150))
ys = map(f->f(x), fs)
finey = map(f->f(finex), fs)

function compare_D_CE_levels(f::Function, x::Array{Float64, 1}, order::Int)

    # Create Dierckx Spline (and data)
    n = length(x)
    y = f(x)
    dspl = Spline1D(x, y; s=0., k=order)

    # Create CompEcon Spline (and data)
    cebasis = Basis(SplineParams(x, 0, order))
    xx = nodes(cebasis)[1]
    yy = f(xx)
    c, bs = funfitxy(cebasis, xx, yy)

    # Evaluate Splines on finer grid
    xfine = collect(linspace(x[1], x[end], 2*n + 1))
    yfine = f(xfine)

    deval = Dierckx.evaluate(dspl, xfine)
    derr = abs(yfine - deval)

    ceeval = CompEcon.funeval(c, cebasis, xfine)
    ceerr = abs(yfine - ceeval)

    # # Evaluate Derivative
    # dderiv1 = Dierckx.derivative(dspl1, finex)
    # cederiv1 = CompEcon.evaluate(ceinterp1, finex, order=1)

    return derr, ceerr
end

function compare_D_CE_deriv(f::Function, fp::Function, x::Array{Float64, 1}, order::Int)

    # Create Dierckx Spline (and data)
    n = length(x)
    y = f(x)
    dspl = Spline1D(x, y; s=0., k=order)

    # Create CompEcon Spline (and data)
    cebasis = Basis(SplineParams(x, 0, order))
    xx = nodes(cebasis)[1]
    yy = f(xx)
    c, bs = funfitxy(cebasis, xx, yy)

    # Evaluate Splines on finer grid
    xfine = collect(linspace(x[1], x[end], 2*n + 1))
    ypfine = fp(xfine)

    devalderiv = Dierckx.derivative(dspl, xfine)
    derr = abs(ypfine - devalderiv)

    ceevalderiv = CompEcon.funeval(c, cebasis, xfine, 1)
    ceerr = abs(ypfine - ceevalderiv)

    return derr, ceerr
end


for i=1:nf

    curr_f, curr_fp = fs[i], fps[i]

    derr1_levels, ceerr1_levels = compare_D_CE_levels(curr_f, x, 1)
    derr2_levels, ceerr2_levels = compare_D_CE_levels(curr_f, x, 2)
    derr3_levels, ceerr3_levels = compare_D_CE_levels(curr_f, x, 3)

    # derr1_deriv, ceerr1_deriv = compare_D_CE_deriv(curr_f, curr_fp, x, 1)
    derr2_deriv, ceerr2_deriv = compare_D_CE_deriv(curr_f, curr_fp, x, 2)
    derr3_deriv, ceerr3_deriv = compare_D_CE_deriv(curr_f, curr_fp, x, 3)

    println("Linear Approximation to Function")
    println("Dierckx Level Errors \n\t Max: $(maximum(derr1_levels)) \n\t Min: $(minimum(derr1_levels)) \n\t Mean: $(mean(derr1_levels))")
    # println("Dierckx Deriv Errors \n\t Max: $(maximum(derr1_deriv)) \n\t Min: $(minimum(derr1_deriv)) \n\t Mean: $(mean(derr1_deriv))")
    println("CompEcon Level Errors \n\t Max: $(maximum(ceerr1_levels)) \n\t Min: $(minimum(ceerr1_levels)) \n\t Mean: $(mean(ceerr1_levels))")
    # println("CompEcon Deriv Errors \n\t Max: $(maximum(ceerr1_deriv)) \n\t Min: $(minimum(ceerr1_deriv)) \n\t Mean: $(mean(ceerr1_deriv))")

    println("Quadratic Approximation to Function")
    println("Dierckx Level Errors \n\t Max: $(maximum(derr2_levels)) \n\t Min: $(minimum(derr2_levels)) \n\t Mean: $(mean(derr2_levels))")
    println("Dierckx Deriv Errors \n\t Max: $(maximum(derr2_deriv)) \n\t Min: $(minimum(derr2_deriv)) \n\t Mean: $(mean(derr2_deriv))")
    println("CompEcon Level Errors \n\t Max: $(maximum(ceerr2_levels)) \n\t Min: $(minimum(ceerr2_levels)) \n\t Mean: $(mean(ceerr2_levels))")
    println("CompEcon Deriv Errors \n\t Max: $(maximum(ceerr2_deriv)) \n\t Min: $(minimum(ceerr2_deriv)) \n\t Mean: $(mean(ceerr2_deriv))")

    println("Cubic Approximation to Function")
    println("Dierckx Level Errors \n\t Max: $(maximum(derr3_levels)) \n\t Min: $(minimum(derr3_levels)) \n\t Mean: $(mean(derr3_levels))")
    println("Dierckx Deriv Errors \n\t Max: $(maximum(derr3_deriv)) \n\t Min: $(minimum(derr3_deriv)) \n\t Mean: $(mean(derr3_deriv))")
    println("CompEcon Level Errors \n\t Max: $(maximum(ceerr3_levels)) \n\t Min: $(minimum(ceerr3_levels)) \n\t Mean: $(mean(ceerr3_levels))")
    println("CompEcon Deriv Errors \n\t Max: $(maximum(ceerr3_deriv)) \n\t Min: $(minimum(ceerr3_deriv)) \n\t Mean: $(mean(ceerr3_deriv))")


end