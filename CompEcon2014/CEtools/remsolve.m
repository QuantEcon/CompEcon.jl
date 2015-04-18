%% REMSOLVE  
%
%  Discrete-Time Continuous-State Continuous-Response Infinite-Horizon
%  Rational Expectations Equilibrium Solver
%
%  Uses the method of collocation and function iteration to compute an
%  approximate solution to the discrete-time infinite-horizon rational
%  expectations model with continuous state s and continuous response x.
%  Returns response function basis function coefficients c and optimal
%  responses x on a refined grid of continuous state nodes.
%
% Usage
%   [c,sr,xr,fr,resid] = remsolve(model,basis,x)
% Let
%   nc = number of state collocation nodes
%   ns = number of state nodes on refined output grid
%   ds = dimension of state s
%   dx = dimension of response x
%   de = dimension of shock e
% Input
%   model : model structure
%   basis : approximation structure
%   x     : nc.dx  initial guess for responses at collocation nodes
% Output
%   c     : nc.1   response function approximant basis function coefficients
%   sr    : ns.ds  refined state nodes
%   xr    : ns.dx  responses at refined state nodes
%   fr    : ns.dx  arbitrage profits at refined state nodes
%   resid : ns.1   response residuals at refined state nodes
%   Note: The number of refined state nodes ns is determined by the grid 
%   refinement factor m, an optional parameter that may be set by the user.  
%   Its default value of 10.  If m>0, the routine will return responses and 
%   residuals on a refined grid of ns uniformly-spaced state nodes that 
%   possesses, along each dimension, m times the coordinates possessed by 
%   the state collocation node grid. If m=0 is set by the user, the routine
%   will return the responses and residuals on the original state 
%   collocation nodes.  
% Model Structure Fields
%   func      : name of function file (see below)
%   params    : parameters passed to function file
%   ds        : dimension of state s
%   dx        : dimension of response x
%   e         : ne.de shocks
%   w         : ne.1  shock probabilities 
% Function File
%   A user-supplied function that returns the bound, arbitrage, and state 
%   transitions, and and their first derivatives with respect to the 
%   response x, at an arbitrary number ns of states s and responses x
%   according to the format
%     [out1,out2,out3,out4] = func(flag,s,x,sn,xn,e,params)
%   Function File Input
%     flag      : flag indicating function to be evaluated 
%     s         : ns.ds states
%     x         : ns.dx responses
%     sn        : ns.ds states next period
%     xn        : ns.dx responses next period
%     e         : ns.de shocks
%     params    : parameters passed to function file
%   Function File Output
%     if flag = 'b', returns lower and upper bounds on response x
%       out1    : ns.dx lower bounds on response x
%       out2    : ns.dx upper bounds on response x
%       out3    : empty
%       out4    : empty
%     if flag = 'f', returns marginal arbitrage profit and derivatives
%       out1    : ns.dx    marginal arbitrage profit
%       out2    : ns.dx.dx first derivative of f with respect to x
%       out3    : ns.dx.ds first derivative of f with respect to sn
%       out4    : ns.dx.dx first derivative of f with respect to xn
%     if flag = 'g', returns state transition and derivatives
%       out1    : ns.ds    state next period g
%       out2    : ns.ds.dx first derivative of g with respect to x
%       out3    : empty
%       out4    : empty
% Options
%   These options may be set by user with OPTSET (defaults in parentheses):
%   tol       : solution algorithm convergence tolerance (sqrt(eps))
%   maxit     : solution algorithm maximum number of iterations (500)
%   nres      : state grid refinement factor (10)
%   limcheck  : 1 if extrapolation to be checked, 0 otherwise

%  Copyright(c) 1997-2014
%    Mario J. Miranda - miranda.4@osu.edu
%    Paul L. Fackler  - paul_fackler@ncsu.edu


function [c,sr,xr,fr,resid] = remsolve(model,basis,x)

% Set options to defaults, if not set by user with OPTSET (see above)
tol       = optget('remsolve','tol',sqrt(eps));
maxit     = optget('remsolve','maxit',500);
nres      = optget('remsolve','nres',10);
limcheck  = optget('remsolve','limcheck',1);

% Set model structure parameters to defaults, if nonexistent
if ~isfield(model,'ds'), model.ds = 1; end % dimension of state s
if ~isfield(model,'dx'), model.dx = 1; end % dimension of response x
if ~isfield(model,'e'),  model.e  = 0; end % shocks
if ~isfield(model,'w'),  model.w  = 1; end % shock probabilities
  
% Unpack model structure
func   = model.func;        % name of function file
params = model.params;      % parameters passed to function file
ds = model.ds;              % dimension of state s
dx = model.dx;              % dimension of response x
e  = model.e;              	% shocks

% Determine number of collocation nodes
n  = basis.n;               % number of state collocation node coordinates by dimension
nc = prod(n);               % number of state collocation nodes

% Set initial guess for values and responses if not provided
if nargin<3, x = zeros(nc,dx); end

% Compute collocation nodes and matrix and basis coefficients
s = funnode(basis);         % state collocation coordinates
Phi = funbasex(basis);     	% collocation matrix
c = funfitxy(basis,Phi,x);	% initial basis coefficients

% Solve arbitrage equlibrium conditions
tic
display('Solving collocation equation by function iteration.')
for it=1:maxit                                   
  cold = c;                                       % store old basis coefficients
  [f,x] = arbit(model,basis,s,x,c);               % update expectation and response variables
  c = funfitxy(basis,Phi,x);                      % update basis coefficient
  change = norm(c-cold,inf);                      % compute change
  fprintf ('%4i %10.1e\n',it,change)              % print change
  if change<tol, break, end;                      % convergence check
end
if it==maxit, display('Failure to converge in remsolve.'), end
fprintf('Elapsed Time = %7.2f Seconds\n',toc)

% Check whether state transitions remain within approximation interval
if limcheck
  snmin =  inf;
  snmax = -inf;
  for k=1:size(e,1);
    kk = k+zeros(nc,1);
    sn = feval(func,'g',s,x,[],[],e(kk,:),params{:});
    snmin = min(snmin,min(sn));
    snmax = max(snmax,max(sn));
  end
  if any(snmin<basis.a-eps), display('Extrapolating below smin. Decrease smin.'), end;
  if any(snmax>basis.b+eps), display('Extrapolating above smax. Increase smax.'), end;
  disp([basis.a' snmin' snmax' basis.b'])
end

% Compute values and responses on refined state grid, if requested.
if nres>0
  nr = nres*n;
  for is=1:ds
    srcoord{is} = linspace(basis.a(is),basis.b(is),nr(is))';
  end
  sr = gridmake(srcoord);
  xr = funeval(funfitxy(basis,Phi,x),basis,sr);   	% rough guess for responses at evaluation points
  [fr,xr] = arbit(model,basis,sr,xr,c);             % shadow prices and responses at evaluation points
else
  sr = s;
  fr = f;
  xr = x;
end

% Compute response residual on output state grid
resid = funeval(c,basis,sr) - xr;


%% ARBIT
%
%   Solves artibrage conditions at ns continuous state nodes s as a
%   nonlinear complementarity problem.  Uses successive linearization, 
%   which requires, at each iteration, evaluation of arbitrage function and 
%   its derivatives with respect to response x.

function [F,x] = arbit(model,basis,s,x,c)

% Set options to defaults, if not set by user with OPTSET (see above)
maxit     = optget('arbit','maxit',300);
tol       = optget('arbit','tol',sqrt(eps));
lcpmethod = optget('arbit','lcpmethod','minmax');

% Unpack model structure
func   = model.func;        % name of function file
params = model.params;      % parameters passed to function file
ds     = model.ds;          % dimension of state s
dx     = model.dx;          % dimension of response x
e      = model.e;           % shocks
w      = model.w;           % shock probabilties

% Solve arbitrage  conditions
[xl,xu] = feval(func,'b',s,x,[],[],[],params{:});
for it=1:maxit
  ns = size(s,1);
  F  = zeros(ns,dx);
  Fx = zeros(ns,dx,dx);
  for k=1:length(w)
    kk = k+zeros(ns,1);
    [sn,snx]  = feval(func,'g',s,x,[],[],e(kk,:),params{:});
    [xn,xnsn] = fund(c,basis,sn,1);
    [f,fx,fsn,fxn] = feval(func,'f',s,x,sn,xn,[],params{:});
    F = F + w(k)*f;
    for j=1:dx
      for ix=1:dx
        for isn=1:ds
          for ixn=1:dx
            Fx(:,j,ix) = Fx(:,j,ix) + w(k)*(fx(:,j,ix)+(fsn(:,j,isn)+fxn(:,j,ixn).*xnsn(:,ixn,isn)).*snx(:,isn,ix));
          end
        end
      end
    end
  end
  [lcpf,delx] = lcpstep(lcpmethod,x,xl,xu,F,Fx);
  x = x + delx;
  if norm(delx(:))< tol, break, end;
end