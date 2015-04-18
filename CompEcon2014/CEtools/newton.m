%% NEWTON   
%
%  Computes root of function via Newton's Method with backstepping
%
%  Usage
%   [x,fval] = newton(f,x,varargin)
%  Input
%    f        : name of function of form [fval,fjac]=f(x,vargin)
%    x        : initial guess for root
%    vargin   : additional arguments for f (optional)
%  Output
%    x        : root of f
%    fval     : function value estimate
%  Options
%    Options may be set by user with OPTSET (defaults in parentheses):
%    tol       : convergence tolerance
%    maxit     : maximum number of iterations
%    maxsteps  : maximum number of backsteps
%    showiters : display results of each iteration

%  Copyright(c) 1997-2010
%   Mario J. Miranda - miranda.4@osu.edu
%   Paul L. Fackler  - paul_fackler@ncsu.edu

function [x,fval] = newton(f,x,varargin)

maxit     = optget('newton','maxit',100);
tol       = optget('newton','tol',sqrt(eps));
maxsteps  = optget('newton','maxsteps',25);
showiters = optget('newton','showiters',0);

for it=1:maxit
   [fval,fjac] = feval(f,x,varargin{:});
   fnorm = norm(fval,inf);
   if fnorm<tol, return, end
   dx = real(-(fjac\fval));
   fnormold = inf;  
   for backstep=1:maxsteps 
      fvalnew = feval(f,x+dx,varargin{:});
      fnormnew = norm(fvalnew,inf);
      if fnormnew<fnorm, break, end
      if fnormold<fnormnew, dx=2*dx; break, end
      fnormold = fnormnew;
      dx = dx/2;          
   end
   x = x+dx;
   if showiters, fprintf('%4i %4i %6.2e\n',[it backstep fnormnew]); end
end
warning('Failure to converge in newton');