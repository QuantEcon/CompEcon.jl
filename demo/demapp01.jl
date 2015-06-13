import CompEcon
using PyPlot
using PyCall

reload("CompEcon")
## DEMAPP01 Approximating functions on R

# This m-file illustrates how to use CompEcon Toolbox routines to construct
# and operate with an approximant for a function defined on an interval of
# the real line.

# In particular, we construct an approximant for f(x)=exp(-x) on the
# interval [-1,1].  The function used in this illustration posseses a
# closed-form, which will allow us to measure approximation error precisely.
# Of course, in practical applications, the function to be approximated
# will not possess a known closed-form.

# In order to carry out the exercise, one must first code the function to
# be approximated at arbitrary points.  The required code is presented at
# the end of this m-file (see below). Let's begin:

# Preliminary tasks
function main()
        f(x) = exp(-x)

    # Set the endpoints of approximation interval:
    a =  -1                            # left endpoint
    b =   1                            # right endpoint

    # Choose an approximation scheme. In this case, let us use an order 10
    # Chebychev approximation scheme:
    n = 10                             # order of approximation
    basis = CompEcon.fundefn(:cheb, n, a, b)      # define basis

    # Compute the basis coefficients c.  There are various way to do this:
    # One may use funfitf:
    c = CompEcon.funfitf(basis, f)

    # ... or one may compute the standard approximation nodes x and
    # corresponding function values y and use funfitxy:
    x = CompEcon.funnode(basis)[1]
    y = f(x)
    c = CompEcon.funfitxy(basis, x, y)

    # ... or one compute the standard approximation nodes x, corresponding
    # function values y, and the interpolation matrix phi, and solve the
    # interpolation equation directly using the backslash operator:
    x = CompEcon.funnode(basis)[1]
    y = f(x)
    phi = CompEcon.funbase(basis)
    c = phi\y

    # Having computed the basis coefficients, one may now evaluate the
    # approximant at any point x using funeval:
    x = [0.0]
    y = CompEcon.funeval(c, basis, x)[1]
    println("The approximate value of exp(-x) at x=0 is $y")
    println("The ''exact'' value of exp(-x) at x=0 is 1")

    # ... one may also evaluate the approximant's first and second derivatives
    # at x:
    d1 = CompEcon.funeval(c, basis, x, 1)[1][1]
    d2 = CompEcon.funeval(c, basis, x, 2)[1][1]
    println("The approximate first derivative of exp(-x) at x=0 is $d1")
    println("The ''exact'' first derivative of exp(-x) at x=0 is -1")
    println("The approximate second derivative of exp(-x) at x=0 is $d2")
    println("The ''exact'' second derivative of exp(-x) at x=0 is 1")

    # ... and one may even evaluate the approximant's definite integral between
    # the left endpoint a and x:
    int = CompEcon.funeval(c, basis, x, -1)[1]
    println("The approximate integral of exp(-x) between x=-1 and x=0 is $int")
    println("The `exact` integral of exp(-x) between x ∈[-1, 0] is $(exp(1)-1)")

    # One may evaluate the accuracy of the Chebychev polynomial approximant by
    # computing the approximation error on a highly refined grid of points:
    ngrid = 5001                       # number of grid nodes
    xgrid = linspace(a, b, ngrid)        # generate refined grid for plotting


    function plot_approx(f::Function, basistype::Symbol, n, a, b, xgrid, k=0)
        basis = basistype == :spli ? CompEcon.fundefn(:spli, n, a, b, k) :
                basistype == :cheb ? CompEcon.fundefn(:cheb, n, a, b) :
                error("Not doing this")

        # compute approximation
        c = CompEcon.funfitf(basis, f)
        yapp = CompEcon.funeval(c, basis, xgrid)[1]

        # compute exact
        yact = f(xgrid)

        # plot error
        fig, ax = plt.subplots()
        ax[:plot](xgrid,yapp-yact)
        ax[:plot](xgrid, zeros(ngrid), "k--", linewidth=2)
        ax[:set_xlabel]("x")
        ax[:set_ylabel]("Error")
        nm = basistype == :spli ? "Degree $k spline" : "Degree $n Chebyshev"
        ax[:set_title]("$nm Approximation Error for exp(-x)")
    end

    plot_approx(f, :cheb, n, a, b, xgrid)
    plot_approx(f, :spli, 21, a, b, xgrid, 3)
    # plot_approx(f, :cheb, 31, a, b, xgrid, 1)
end

main()
