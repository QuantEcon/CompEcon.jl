import CompEcon
using PyPlot
using PyCall

function main()
    f1(x) = exp(-2x)
    g1(x) = -2exp(-2x)

    function univariate_interp(f, g)
        # fit approximant
        n = 10; a = -1; b = 1
        basis = CompEcon.fundefn(:cheb, n, a, b)
        c = CompEcon.funfitf(basis, f1)

        # graph approximation error for function and derivative
        x = CompEcon.nodeunif(1001, a, b)[1]
        yact = f(x)
        dact = g(x)
        yfit = CompEcon.funeval(c, basis, x)[1]
        dfit = CompEcon.funeval(c, basis, x, 1)[1]


        # Nice plot of function approximation error
        fig, ax = plt.subplots()
        ax[:plot](x, yfit-yact)
        ax[:plot](x, 0*x, "k--", linewidth=2)
        ax[:set_xlabel]("x")
        ax[:set_ylabel]("Error")
        ax[:set_title](L"Error approximating function using Chebyshev")

        # Nice plot of derivative approximation error
        fig, ax = plt.subplots()
        ax[:plot](x, dfit-dact);
        ax[:plot](x, 0*x, "k--", linewidth=2)
        ax[:set_xlabel]("x")
        ax[:set_ylabel]("Error")
        ax[:set_title](L"Error approximating derivative using Chebyshev")
    end


    f2(x) = cos(x[:, 1]) ./ exp(x[:, 2])
    function bivaraite_interp(f)
        n = [7 7]
        a = [0 0]
        b = [1 1]
        basis = CompEcon.fundefn(:cheb, n, a, b)
        c = CompEcon.funfitf(basis, f2)

        fig = plt.figure()
        nplot = [101 101]
        x, xcoord = CompEcon.nodeunif(nplot, a, b)
        yfit = CompEcon.funeval(c, basis, x)[1]
        error = reshape(yfit - f2(x), nplot...)
        poly3d = surf(xcoord[1], xcoord[2], error)
        ax = poly3d[:get_axes]()
        ax[:set_title]("Chebyshev approximation error")

    end
end

main()
