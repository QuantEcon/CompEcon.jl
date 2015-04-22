# CompEcon

[![Build Status](https://travis-ci.org/spencerlyon2/CompEcon.jl.svg?branch=master)](https://travis-ci.org/spencerlyon2/CompEcon.jl)



## Quick (and incomplete intro)

Right now you can use two versions of these functions:

- original Matlab API from Miranda and Fackler
- Type based Julia API

### Matlab-esque interface

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
c = funfitxy(basis, x, y)

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

### Julian API

The corresponding code with the more Julian API currently (2015-04-22 17:45:10) looks like this:

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
basis = Basis(Cheb, n, a, b)     # define basis

# Compute the basis coefficients c.  There are various way to do this:
# One may use funfitf:
c = funfitf(basis, f)

# ... or one may compute the standard approximation nodes x and corresponding
# function values y and use funfitxy:
x = nodes(basis)[1]
y = f(x)
c = funfitxy(basis, x, y)[1]

# ... or one compute the standard approximation nodes x, corresponding
# function values y, and the interpolation matrix phi, and solve the
# interpolation equation directly using the backslash operator:
x = nodes(basis)[1]
y = f(x)
phi = BasisStructure(basis).vals[1]
c = phi\y

# Having computed the basis coefficients, one may now evaluate the
# approximant at any point x using funeval:
x = [0.0]
y = funeval(c, basis, x)[1]
```


The Julia style API version doesn't look any simpler than the original API, but users are still strongly encouraged to use it. The reason for this is that the original API code will not change much, if any, in the future. However, the Julia-style code will be updated and improved.

## Example

Here's another example of how to use the Julia-based API to set up multi-dimensional basis structure and do some things with it.

```julia
ygrid0 = linspace(-4, 4, 10)
agrid0 = linspace(0.0.^0.4, 100.0.^0.4, 25).^(1/0.4)

# method one, using the Basis constructor multiple times
basis = Basis(Basis(Spline, agrid0, 0, 3),  # cubic spline
              Basis(Spline, ygrid0, 0, 1))  # linear

# method two, constructing separately, then calling `Basis` with the two
a_basis = Basis(Spline, agrid0, 0, 3)
y_basis = Basis(Spline, ygrid0, 0, 1)
basis = Basis(a_basis, y_basis)

# could also use unicode `\times[TAB]` to combine basis
basis = a_basis × y_basis

# Construct state vector (matrix). Note that splidef (called by
# fundef) adds breakpoints to the original grid we gave it, so let's
# extract the actual grid points used from the second argument
S, (agrid, ygrid) = nodes(basis)

# construct basis matrix and its lu-factorization for very fast inversion
# NOTE: I am doing this in a round-about way. I could have just done
#       Φ = BasisStructure(basis), but doing it this way gives me the direct
#       representation so I get Φ_y without repeating any calculations
Φ_direct = BasisStructure(basis, S, [0 0], Direct())
Φ_y = Φ_direct.vals[2]
Φ = convert(Expanded, Φ_direct, [0 0]).vals[1]
lu_Φ = lufact(Φ)
```

## Notes

Tests still need to be written.

We will probably add another type `Interpoland` that carries around at least a `Basis` and the coefficient vector to make it easier to call `funeval` without having to keep track coefficients separately.
