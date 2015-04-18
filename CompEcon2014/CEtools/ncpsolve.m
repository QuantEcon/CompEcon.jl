%% NCPSOLVE  
%
%  Solves the general d-dimensional nonlinear complementatity problem
%     a_i <= x_i <= b_i, i=1,2,...,d
%     x_i > a_i => f_i(x) => 0, i=1,2,...,d
%     x_i < b_i => f_i(x) =< 0, i=1,2,...,d
%  where f is a function from R^d to R^d and a and b are d by 1 vectors 
%  with a<=b. Problem is solved by applying a safeguarded Newton method to
%  the equivalent re-forrmulated minimax or semismooth rootfinding problem.
%
%  Usage
%    [x,fval] = ncpsolve(func,a,b,x,<params>)
%  Input
%    func      : name of user-supplied function (see below)
%    a         : d by 1 vector, left bound on x
%    b         : d by 1 vector, right bound on x
%    x         : d by 1 vector, initial guess for solution
%    params    : optional parameters passed to func
%  Output
%    x         : d by 1 solution to ncp
%    fval      : d by 1 function value at x
%  Function File
%    func.m a user-supplied function that returns the d by 1 vector of
%    values and d by d Jacobian of the function f according to the format
%      [fval,fjac] = func(x,<params>)
%    function input
%      x         : d by 1 vector
%      <params>  : optional parameters passed fo func
%    function output
%      fval      : d by 1 vector of function values
%      fjac      : d by d matris of partial derivatives
%  Options
%    Options may be set by user with OPTSET (defaults in parentheses):
%    tol       : convergence tolerance (sqrt(eps)_
%    maxit     : maximum number of iterations (100)
%    maxsteps  : maximum number of backsteps (10)
%    showiters : display results of each iteration (0)
%    type      : rootproblem transform, 'ssmooth' or 'minmax' (ss)

%  Copyright(c) 1997-2010
%   Mario J. Miranda - miranda.4@osu.edu
%   Paul L. Fackler  - paul_fackler@ncsu.edu

function [x,ftmp,ierr] = ncpsolve(f,a,b,x,varargin)

maxit     = optget('ncpsolve','maxit',100);
tol       = optget('ncpsolve','tol',sqrt(eps));
maxsteps  = optget('ncpsolve','maxsteps',10);
showiters = optget('ncpsolve','showiters',0);
type      = optget('ncpsolve','type','ssmooth');

if nargin < 4, x=zeros(size(a)); end

for it=1:maxit
  [fval,fjac] = feval(f,x,varargin{:});
  if isempty(fjac)
     fjac = fdjac(f,x,varargin{:});
  end
  [ftmp,fjac] = feval(type,x,a,b,fval,fjac);
  fnorm = norm(ftmp,inf);
  if fnorm<tol, break, end
  dx = -(fjac\ftmp);
  fnormold = inf;
  for backstep=1:maxsteps
    xnew = x + dx;
    fnew = feval(f,xnew,varargin{:});
    fnew = feval(type,xnew,a,b,fnew);
    fnormnew = norm(fnew,inf);
    if fnormnew<fnorm, break, end
    if fnormold<fnormnew, dx=2*dx; break, end
    fnormold = fnormnew;
    dx = dx/2;
  end
  x = x+dx;
  if showiters, fprintf('%4i %4i %6.2e\n',[it backstep fnormnew]); end
end
if it==maxit
  ierr = 1;
  if nargout<3
    warning('Failure to converge in ncpsolve')
  end
else
  ierr = 0;
end
x = real(x);
x = max(a,x);
x = min(b,x);