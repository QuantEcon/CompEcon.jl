using CompEcon
using Dierckx


# Create functions to test against
f1(x) = sin(x)
f2(x) = exp(x)
f3(x) = 1.5 + log(x)
fs = [f1, f2, f3]
nf = length(fs)

x = collect(linspace(1e-4, 2*pi, 35))
finex = collect(linspace(1e-2, 2*pi - 1e-2, 150))
ys = map(f->f(x), fs)
finey = map(f->f(finex), fs)

for i=1:nf

    #
    # Check Linear
    #
    dspl1 = Spline1D(x, ys[i]; s=0., k=1)
    cebasis1 = Basis(SplineParams(x, 0, 1))
    c1, bs1 = funfitxy(cebasis1, x, ys[i])

    # Evaluate Function
    deval1 = Dierckx.evaluate(dspl1, finex)
    derr1 = abs(finey[i] - deval1)
    ceeval1 = CompEcon.funeval(c1, cebasis1, finex)
    ceerr1 = abs(finey[i] - ceeval1)

    # Evaluate Derivative
    dderiv1 = Dierckx.derivative(dspl1, finex)
    cederiv1 = CompEcon.evaluate(ceinterp1, finex, order=1)


    println("Linear Approximation to Function")
    println("--------------------")
    println("Dierckx Errors \n\t Max: $(maximum(derr1)) \n\t Min: $(minimum(derr1)) \n\t Mean: $(mean(derr1))")
    println("CompEcon Errors \n\t Max: $(maximum(ceerr1)) \n\t Min: $(minimum(ceerr1)) \n\t Mean: $(mean(ceerr1))")

    #
    # Check Quadratic
    #
    dspl2 = Spline1D(x, ys[i]; s=0., k=2)
    cebasis2 = Basis(SplineParams(x, 0, 2))
    c2, bs2 = funfitxy(cebasis2, x, ys[i])

    # Evaluate Function
    deval2 = Dierckx.evaluate(dspl2, finex)
    derr2 = abs(finey[i] - deval1)
    ceeval2 = CompEcon.funeval(c2, cebasis2, finex)
    ceerr2 = abs(finey[i] - ceeval2)

    # Evaluate Derivative
    dderiv2 = Dierckx.derivative(dspl2, finex)
    cederiv2 = CompEcon.evaluate(ceinterp2, finex, order=1)


    println("Quadratic Approximation to Function")
    println("--------------------")
    println("Dierckx Errors \n\t Max: $(maximum(derr2)) \n\t Min: $(minimum(derr2)) \n\t Mean: $(mean(derr2))")
    println("CompEcon Errors \n\t Max: $(maximum(ceerr2)) \n\t Min: $(minimum(ceerr2)) \n\t Mean: $(mean(ceerr2))")

    #
    # Check Cubic
    #
    dspl3 = Spline1D(x, ys[i]; s=0., k=3)
    cebasis3 = Basis(SplineParams(x, 0, 3))
    c3, bs3 = funfitxy(cebasis3, x, ys[i])


    # Evaluate
    deval3 = Dierckx.evaluate(dspl3, finex)
    derr3 = abs(finey[i] - deval3)
    ceeval3 = CompEcon.funeval(c3, cebasis3, finex)
    ceerr3 = abs(finey[i] - ceeval3)
    println("Cubic Approximation")
    println("--------------------")
    println("Dierckx Errors \n\t Max: $(maximum(derr3)) \n\t Min: $(minimum(derr3)) \n\t Mean: $(mean(derr3))")
    println("CompEcon Errors \n\t Max: $(maximum(ceerr3)) \n\t Min: $(minimum(ceerr3)) \n\t Mean: $(mean(ceerr3))")
end