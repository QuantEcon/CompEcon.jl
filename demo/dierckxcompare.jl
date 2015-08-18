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
    println("Linear Approximation")
    println("--------------------")
    dspl1 = Spline1D(x, ys[i]; s-0., k=1)
    cebasis1 = Basis(Spline, x, 0, 1)
    cecoeffs1, cebs1 = funfitxy(Basis(Spline, x, 0, 1), x, ys[i])
    ceinterp1 = CompEcon.Interpoland(cebasis1, cecoeffs1, cebs1)

    # Evaluate
    deval3 = Dierckx.evaluate(dspl3, finex)
    derr3 = abs(finey[i] - deval3)
    ceeval3 = CompEcon.evaluate(ceinterp3, finex)
    ceerr3 = abs(finey[i] - ceeval3)
    println("Dierckx Errors \n\t Max: $(maximum(derr1)) \n\t Min: $(minimum(derr1)) \n\t Mean: $(mean(derr1))")
    println("CompEcon Errors \n\t Max: $(maximum(ceerr1)) \n\t Min: $(minimum(ceerr1)) \n\t Mean: $(mean(ceerr1))")

    #
    # Check Cubic
    #
    println("Cubic Approximation")
    println("--------------------")
    dspl3 = Spline1D(x, ys[i]; s=0., k=3)
    cebasis3 = Basis(Spline, x, 0, 3)
    cecoeffs3, cebs3 = funfitxy(cebasis3, x, ys[i])
    ceinterp3 = CompEcon.Interpoland(cebasis3, cecoeffs3, cebs3)

    # Evaluate
    deval3 = Dierckx.evaluate(dspl3, finex)
    derr3 = abs(finey[i] - deval3)
    ceeval3 = CompEcon.evaluate(ceinterp3, finex)
    ceerr3 = abs(finey[i] - ceeval3)
    println("Dierckx Errors \n\t Max: $(maximum(derr3)) \n\t Min: $(minimum(derr3)) \n\t Mean: $(mean(derr3))")
    println("CompEcon Errors \n\t Max: $(maximum(ceerr3)) \n\t Min: $(minimum(ceerr3)) \n\t Mean: $(mean(ceerr3))")