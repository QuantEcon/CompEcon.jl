[![Build Status](https://travis-ci.org/QuantEcon/CompEcon.jl.svg?branch=master)](https://travis-ci.org/QuantEcon/CompEcon.jl) [![codecov.io](http://codecov.io/github/QuantEcon/CompEcon.jl/coverage.svg?branch=master)](http://codecov.io/github/QuantEcon/CompEcon.jl?branch=master) [![Coverage Status](https://coveralls.io/repos/QuantEcon/CompEcon.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/QuantEcon/CompEcon.jl?branch=master)

This package is a Julia implementation of the routines originally contained in the [CompEcon Matlab toolbox](http://www4.ncsu.edu/~pfackler/compecon/toolbox.html) by Paul Fackler and Mario Miranda. The original Matlab code was written to accompany the publication

> Miranda, Mario J., and Paul L. Fackler. Applied computational economics and finance. MIT press, 2004.

This work is derivative of their work and has been licensed with their permission.

# CompEcon

This package is a wrapper around [BasisMatrices.jl](https://github.com/QuantEcon/BasisMatrices.jl) and provides an API similar to the original [CompEcon matlab library](http://www4.ncsu.edu/~pfackler/compecon/toolbox.html) by Miranda and Fackler. For best use of the underlying routines, we recommend using the BasisMatrices.jl API.

The Matlab style API here is as close to the original library as possible (differences are based mostly on syntax). To see what this means, consider the following Matlab example (taken from `demapp01.m`):

```matlab
% function to approximate
f = @(x) exp(-x)

% Set the endpoints of approximation interval:
a =  -1;                            % left endpoint
b =   1;                            % right endpoint

% Choose an approximation scheme. In this case, let us use an order 10
% Chebychev approximation scheme:
n = 10;                             % order of approximation
basis = fundefn('cheb',n,a,b);      % define basis

% Compute the basis coefficients c.  There are various way to do this:
% One may use funfitf:
c = funfitf(basis,@f);

% ... or one may compute the standard approximation nodes x and corresponding
% function values y and use funfitxy:
x = funnode(basis);
y = f(x);
c = funfitxy(basis,x,y);

% ... or one compute the standard approximation nodes x, corresponding
% function values y, and the interpolation matrix phi, and solve the
% interpolation equation directly using the backslash operator:
x = funnode(basis);
y = f(x);
phi = funbase(basis);
c = phi\y;

% Having computed the basis coefficients, one may now evaluate the
% approximant at any point x using funeval:
x = 0;
y = funeval(c,basis,x);
```


The corresponding Julia code is

```julia
using CompEcon
# function to approximate
f(x) = exp(-x)

# Set the endpoints of approximation interval:
a =  -1                            # left endpoint
b =   1                            # right endpoint

# Choose an approximation scheme. In this case, let us use an order 10
# Chebychev approximation scheme:
n = 10                             # order of approximation
basis = fundefn(:cheb, n, a, b)      # define basis

# Compute the basis coefficients c.  There are various way to do this:
# One may use funfitf:
c = funfitf(basis, f)

# ... or one may compute the standard approximation nodes x and corresponding
# function values y and use funfitxy:
x = funnode(basis)[1]
y = f(x)
c = funfitxy(basis, x, y)[1]

# ... or one compute the standard approximation nodes x, corresponding
# function values y, and the interpolation matrix phi, and solve the
# interpolation equation directly using the backslash operator:
x = funnode(basis)[1]
y = f(x)
phi = funbase(basis)
c = phi\y

# Having computed the basis coefficients, one may now evaluate the
# approximant at any point x using funeval:
x = [0.0]
y = funeval(c, basis, x)[1]
```


The main differences are:

- The Julia code uses symbols instead of strings to specify basis functions and refer to objects in the basis structure. The Matlab uses string (we see this above with use of `'cheb'` in Matlab and `:cheb` in Julia)
- The Matlab code relies heavily on the use of `varargout` to only return some objects. The Julia code always returns all objects the Matlab ones might ever return, so we need to be careful about keeping only some of the return arguments. (notice in the calls to `funnode`  and `funeval` we just keep the first output in Julia).
