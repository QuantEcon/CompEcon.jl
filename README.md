# CompEcon

[![Build Status](https://travis-ci.org/spencerlyon2/CompEcon.jl.svg?branch=master)](https://travis-ci.org/spencerlyon2/CompEcon.jl) [![codecov.io](http://codecov.io/github/spencerlyon2/CompEcon.jl/coverage.svg?branch=master)](http://codecov.io/github/spencerlyon2/CompEcon.jl?branch=master) [![Coverage Status](https://coveralls.io/repos/spencerlyon2/CompEcon.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/spencerlyon2/CompEcon.jl?branch=master)



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

The Julia style API version doesn't look any simpler than the original API, but users are still strongly encouraged to use it. The reason for this is that the original API code will not be improved much, if any, in the future. However, we anticipate making significant improvements to the Julian api.

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

## Basic Overview of Julian API

This section provides a sketch of the type based Julian API.

### Theoretical Foundation

To understand the Julian API and type system, we first need to understand the fundamental theory behind the interpolation scheme implemented here. Interpolation in CompEcon is built around three key concepts:

1. An functional `Basis`: for each dimension, the basis specifies
    - family of basis function (B spline, Chebyshev polynomials, ect.)
    - domain (bounds)
    - interpolation nodes (grid on domain)
2. A `BasisStructure`:
    - Represents the evaluation of basis functions at the interpolation nodes
    - Constructed one dimension at a time, then combined with tensor product
3. A coefficient vector: used to map from domain of the `Basis` into real line

### Core types

Functionality implemented around 5 core types (or type families) that relate closely to the theoretical concepts outlined above.

#### Representing the `Basis`

The first two groups of type are helper types used to facilitate construction of the `Basis`. They are the `BasisFamily` and the `BasisParams` types:

First is the `BasisFamily`:

```julia
abstract BasisFamily
immutable Cheb <: BasisFamily end
immutable Lin <: BasisFamily end
immutable Spline <: BasisFamily end

abstract BasisParams
type ChebParams <: BasisParams
    n::Int
    a::Float64
    b::Float64
end

type SplineParams <: BasisParams
    breaks::Vector{Float64}
    evennum::Int
    k::Int
end

type LinParams <: BasisParams
    breaks::Vector{Float64}
    evennum::Int
end
```

`BasisFamily` is an abstract type, whose subtypes are singletons that specify the class of functions in the basis.

`BasisParams` is an abstract type, whose subtypes are type types that hold all information needed to construct the Basis of a particular class

Then we have the central `Basis` type:

```julia
type Basis{N}
    basistype::Vector{BasisFamily}  # Basis family
    n::Vector{Int}                  # number of points and/or basis functions
    a::Vector{Float64}              # lower bound of domain
    b::Vector{Float64}              # upper bound of domain
    params::Vector{BasisParams}     # params to construct basis
end
```

Each field in this object is a vector. The `i`th element of each vector is the value that specifies the commented description for the `i`th dimension.

The `Basis` has support for the following methods:

- A whole slew of constructors
- Conversion to and from the `Dict`s used in the old Matlab-eqsue API
- `getindex(b::Basis, i::Int)`: which extracts the univariate `Basis` along the `i`th dimension
- `ndims`: The number of dimensions
- `length`: the product of the `n` field
- `size(b::Basis, i::Int)`: The `i`th element of the `n` field (number of basis functions in dimension `i`)
- `size(b::Basis)`: `b.n` as a tuple instead of a vector (similar to `size(a::Array)`)
- `==`: test two basis for equality
- `nodes(b::Basis)->(Matrix, Vector{Vector{Float64}})`: the interpolation nodes. the first element is the tensor product of all dimensions, second element is a vector of vectors, where the `i`th element contains the nodes along dimension `i`.

#### `BasisStructure` representation

Next we turn to representing the `BasisStructure`, which is responsible for keeping track of the basis functions evaluated at the interpolation nodes. To keep traack of this representation, we have another family of helper types:

```julia
abstract AbstractBasisStructureRep
typealias ABSR AbstractBasisStructureRep

immutable Tensor <: ABSR end
immutable Direct <: ABSR end
immutable Expanded <: ABSR end
```

`AbstractBasisStructureRep` is an abstract types, whose subtypes are singletons that specify how the basis matrices are stored. To understand how they are different, we need to see the structure of the `BasisStructure` type:

```julia
type BasisStructure{BST<:ABSR}
    order::Matrix{Int}
    vals::Array{AbstractMatrix}
end
```

The `order` field keeps track of what order of derivative or integral the arrays inside `vals` correspond to.


The content inside `vals` will vary based on the type Parameter `BST<:AbstractBasisStructureRep`:

1. for `BST==Tensor` `vals` will store the evaluation of the basis functions at each of the integration nodes indpendently. Thus, the `[d, i]` element will be the derivative order `d` Basis at the interpolation nodes along the `i`th dimension (each column is a basis function, each row is an interpolation node). This is the most compact and memory efficient representation
2. For `BST==Direct` `vals` will expand along the first dimension (rows) of the array so that for each `i`, the `[d, i]` element will have `length(basis)` rows and `basis.n[i]` (modulo loss or addition of basis functions from derivative/intergral operators.)
3. For `BST==Expanded` `vals` will be expanded along both the rows and the columns and will contain a single `Matrix` for each desired derivative order. This format is the least memory efficient, but simplest conceptually for thinking about how to compute a coefficient vector (if `y` is `f(x)` then `coefs[d] = b.vals[d] \ y`)

#### Convenience `Interpoland` type

Finally the convenient `Interpoland` type:

```julia
type Interpoland{T<:FloatingPoint,N,BST<:ABSR}
    basis::Basis{N}
    coefs::Vector{T}
    bstruct::BasisStructure{BST}
end
```

This type doesn't do a whole lot. It's main purpose is to keep track of the coefficient vector and the `Basis` so the user doesn't have to carry both of them around. It also holds a `BasisStructure` for the evaluation of the basis matrices at the interpolation nodes. This means that if the coefficient vector needs to be updated, this `BasisStructure` will not be re-computed.
