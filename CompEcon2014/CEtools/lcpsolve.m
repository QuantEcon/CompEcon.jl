%% LCPSOLVE 
%
%  Uses safeguarded Newton's method to solve linear complementarity problem
%     z = M*x + q
%     a <= x <= b
%     x_i > a_i => z_i => 0
%     x_i < b_i => z_i =< 0
%
%  Usage
%    [x,z] = lcpsolve(M,q,a,b,x)
%  Input
%    M       : n by 1 matrix
%    q       : n by 1 vector
%    a       : n by 1 vector, left bound on x
%    b       : n by 1 vector, right bound on x
%    x       : n by 1 vector, initial guess for solution
%  Output
%    x       : solution to lcp
%    z       : function value at x
%  Options
%    Options may be set with OPTSET
%    tol     : convergence tolerance
%    maxit   : maximum number of iterations
%    maxsteps: maximum number of backsteps
%    type    : rootproblem transform, 'ssmooth' or 'minmax'

%  Copyright(c) 1997-2010
%   Mario J. Miranda - miranda.4@osu.edu
%   Paul L. Fackler  - paul_fackler@ncsu.edu

function [x,z] = lcpsolve(M,q,a,b,x)

maxit    = optget('lcpsolve','maxit',100);
tol      = optget('lcpsolve','tol',sqrt(eps));
maxsteps = optget('lcpsolve','maxsteps',20);
type     = optget('lcpsolve','type','ssmooth');

if nargin < 5, x=q; end

for it=1:maxit
  z = M*x+q;
  [ztmp,Mtil] = feval(type,x,a,b,z,M);
  fnorm = norm(ztmp);
  if fnorm<tol, return, end
  dx = -(Mtil\ztmp);
  fnormold = inf;
  for backstep=1:maxsteps
    xnew = x + dx;
    znew = M*xnew+q;
    ztmp = feval(type,xnew,a,b,znew);
    fnormnew = norm(ztmp);
    if fnormnew<fnorm, break, end
    if fnormold<fnormnew, dx=dx*2; break, end
    fnormold = fnormnew;
    dx = dx/2;
  end
  x = x+dx;
end
warning('In lcpsolve" failure to converge.')
z = M*x+q;