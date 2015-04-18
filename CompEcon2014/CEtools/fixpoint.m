%% FIXPOINT
%
%  Computes fixed-point of function using function iteration
% 
%  Usage
%    x = fixpoint(g,x,varargin)
%  Input
%    g        : name of function of form gval=g(x)
%    x        : initial guess for fixed-point
%    varargin : parameters passed to function g
%  Output
%    x        : fixed-point of g
%    gval     : function value at x
%  Options
%    maxit    : maximum number of iterations (100)
%    tol      : convergence tolerance (sqrt(eps))

%  Copyright(c) 1997-2010
%   Mario J. Miranda - miranda.4@osu.edu
%   Paul L. Fackler  - paul_fackler@ncsu.edu

function x = fixpoint(g,x,varargin)

% Set options to defaults, if not set by user with OPTSET (see above)
maxit = optget('fixpoint','maxit',100);
tol   = optget('fixpoint','tol',sqrt(eps));

for it=1:maxit
  xold = x;
  x = feval(g,x,varargin{:});
  if norm(x-xold)<tol, return, end
end
warning('Failure to converge in fixpoint')