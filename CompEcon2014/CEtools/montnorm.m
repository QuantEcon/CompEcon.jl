% MONTNORM Computes pseudo-random multivariate normal variates
% USAGE
%   x = montnorm(n,mu,var)
% INPUTS
%   n   : number of variates
%   mu  : 1 by d mean vector (default=0)
%   var : d by d positive definite covariance matrix (default=I)
% OUTPUTS
%   x   : n by d matrix of evaluation nodes
% 
% To compute expectation of f(x), where x is N(mu,var), write a
% function f that returns m-vector of values when passed an m by d
% matrix, and write x=montnorm(n,mu,var); Ef=f(x)/n;

% Copyright (c) 1997-2000, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function x = montnorm(n,mu,var)

if (nargin<2) | isempty(mu)
  d=1; 
  mu = 0;
else
  d = length(mu);
  mu=mu(:)';
end
if (nargin<3) | isempty(var), var = speye(d,d); end

z = randn(n,d);
if d==1
   x = mu+sqrt(var)*z;
else
   x = mu(ones(1,prod(n)),:) + z*chol(var);
end
