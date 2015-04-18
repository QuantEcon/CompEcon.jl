%% QNWLOGN
%
%  Generates Gaussian quadrature nodes and probability weights for 
%  multivariate Lognormal distribution with parameters mu and var, or, 
%  equivalently, for a random vector whose natural logarithm is normally 
%  distributed with mean vector mu and variance matrix var.
%
%  Usage
%    [x,w] = qnwlogn(n,mu,var)
%  Input
%    n   : 1.d number of nodes per dimension
%    mu  : 1.d mean vector (default: zeros)
%    var : d.d positive definite variance matrix (default: identity matrix)
%  Output
%    x   : prod(n).d quadrature nodes
%    w   : prod(n).1 probability weights
%  Note
%    The lognormal distribution is defined on (0,inf)^d.  The mean and
%    and variance of the univariate distribution are exp(mu+var/2) and 
%    (exp(var)-1)exp(2mu+var), respectively.
%  Note
%    To compute Ef(X) when f is real-values and X is Lognormal(mu,var) on
%    R^d, write a Matlab function f that returns an m.1 vector when passed 
%    an m.d matrix, and write [x,w]=qnwnorm(n,mu,var); Ef=w'*f(x).
%  Uses
%    qnwnorm

%  Copyright(c) 1997-2010
%   Mario J. Miranda - miranda.4@osu.edu
%   Paul L. Fackler  - paul_fackler@ncsu.edu

function [x,w] = qnwlogn(n,mu,var)

d = length(n);
if nargin<2, mu  = zeros(1,d); end
if nargin<3, var = eye(d); end
if size(mu,1)>1, mu=mu'; end

[x,w] = qnwnorm(n,mu,var);
x = exp(x);