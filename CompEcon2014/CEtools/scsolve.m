%% SCSOLVE
%
%  General Continuous-Time Stochastic Control Solver
%
%  Solves the continuous time Jacoby-Bellman equation
%   rho*V = max_x f(s,x) + V_s(s)g(s,x) + 0.5trace(sigma(s)*sigma(s)'*V_ss(s))
%  using policy function iteration
%
%  Usage
%    [c,scoord,v,x,resid] = scsolve(model,basis,v)
%  Let
%    ds = dimension of state s
%    dx = dimension of action x
%    n  = number of basis functions per state dimension (1.ds)
%    nc = number of nodes on state collocation grid (prod(n))
%    ns = number of state nodes output
%  Input
%    model  : structured array containing model specifications (see below)
%    basis  : ds-dimensional basis for functions defined on state space
%    v      : nc.1 initial guess for values at collocation nodes 
%            (optional, default is array of zeros)
%  Output
%    c      : nc.1 value function approximant basis function coefficients
%    scoord : ns.ds refined state node grid in coordinate form
%    v      : ns.1 values at refined state nodes
%    x      : ns.dx optimal actions at refined state nodes
%    resid  : ns.1 Bellman equation residuals at refined state nodes
%    Note: If the residual is requested, then values, optimal controls, and
%    residuals are returned on a refined grid of ns state nodes.  The
%    degree of refinement is governed by nr, an optional parameter with
%    default value of 10 that may be set by the user using optset.  If
%    nr>0, the refined state node grid is created by forming the Cartesian
%    product of nr*n equally-spaced coordinates along each state dimension.
%    If nr=0, the routine will return the values, optimal controls, and
%    residuals at the original collocation state nodes. On output, the
%    singleton dimensions of c, v, x, and resid are eliminated.
%  Model Structure
%    The structured array "model" contains different fields that specify
%    essential features of the model to be solved (default values for
%    optional fields in parentheses):
%      func    : name of function file (required)
%      params  : model parameters required by function file (empty) 
%      rho     : discount rate (required)
%  Function File
%    User-supplied function that returns the bound, reward, and transition
%    functions and derivatives with respect to the action x at an arbitrary
%    number of ns states and actions according to the format:
%      [out1,out2,out3] = func(flag,s,x,Vs,<params>)
%    Function File Input
%      flag    : flag that specifies function to be evaluated
%      s       : ns.ds state nodes
%      x       : ns.dx actions
%      Vs      : ns.ds first derivatives of value function
%      params  : user-supplied list of function parameters
%    Function File Output
%      flag = 'x' returns 
%        out   : ns.dx optimal controls at nodes
%      flag = 'f' returns
%        out   : ns.1 reward function values at nodes
%      flag = 'g' returns
%        out   : ns.ds state drift function values at nodes
%      flag = 'sigma' returns 
%        out   : ns.ds.ds state diffusion function values at nodes
%  Options
%    maxit     : maximum number of iterations, collocation equation
%                solution algorithm (500)
%    tol       : convergence tolerance (square root of machine epsilon)
%    showiters : 0/1, 1 to display iteration results
%    nr        : state grid refinement factor (10)
%    output    : whether to generate printed output (1)

% Copyright (c) 1997-2014, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function [c,scoord,v,x,resid] = scsolve(model,basis,v)

% Set user options to defaults, if not set by user with OPTSET (see above)
maxit      = optget('scsolve','maxit',500);
tol        = optget('scsolve','tol',sqrt(eps));
nr         = optget('scsolve','nr',10);
showiters  = optget('scsolve','showiters',1);

% Set model fields to default if nonexistent (see above)
if ~isfield(model,'func'),   error('Missing Function File');  end
if ~isfield(model,'params'), error('Missing Parameter List'); end
if ~isfield(model,'rho'),    error('Missing Discount Rate');  end

% Unpack model structure
func    = model.func;
params  = model.params;
rho     = model.rho;

% Initialize state collocation nodes
ds = basis.d;               % dimension of state space
n  = basis.n;               % number of state collocation nodes by dimension
nc = prod(n);               % total number of state collocation nodes
[s,scoord] = funnode(basis);% collocation nodes

% Compute part of collocation matrix that does not depend on x
bases = funbasex(basis,scoord,[0;1;2]);
B = ctbasemake(rho,bases,0);
sigma = feval(func,'sigma',s,[],[],params{:});
if ~isempty(sigma)
  B = B-ctbasemake(sigma,bases,2);
end

% Store Phi1 in expanded form as it is used repeatedly
Phi1 = ctbasemake([],bases,1);

% Initialize coefficient vector
if nargin<3 || isempty(v)
  c = zeros(nc,1);
  v = zeros(nc,1);
else
  c = ckronxi(bases.vals(1,:),v,ds:-1:1);
end

% Policy function iteration loop
for it=1:maxit
  v0 = v;
  Vs = zeros(size(Phi1{1},1),ds);
  for i=1:ds, Vs(:,i) = Phi1{i}*c; end
  x = feval(func,'x',s,[],Vs,params{:});
  f = feval(func,'f',s,x,[],params{:});
  g = feval(func,'g',s,x,[],params{:});
  c = (B-ctbasemake(g,Phi1,1))\f;
  v = funeval(c,basis,bases);
  if any(isnan(v) | isinf(v))
    error('NaNs or Infs encountered');
  end
  e = max(abs(v-v0));
  if showiters, fprintf('%3i %12.4e\n',it,e); end
  if e<tol, break; end
end
if it>=maxit
  disp(['Algorithm did not converge. Maximum error: ' num2str(e)]);
end

if nargout>4 && nr>0
  if ds==1,
    n = nr*n+1;
    scoord = linspace(basis.a,basis.b,n)';
  else
    for i=1:ds
      n(i) = nr*n(i)+1;
      scoord{i} = linspace(basis.a(i),basis.b(i),n(i))';
    end
  end
  bases = funbasex(basis,scoord,[0;1;2]);
  v = funeval(c,basis,bases);
  Vs = squeeze(funeval(c,basis,bases,eye(ds)));
  s = gridmake(scoord);
  x = feval(func,'x',s,[],Vs,params{:});
  if nargout>4
    f = feval(func,'f',s,x,[],params{:});
    resid = rho*v-f;
    g = feval(func,'g',s,x,[],params{:});
    resid = resid-sum(Vs.*g,2);
    sigma = feval(func,'sigma',s,[],[],params{:});
    if ~isempty(sigma)
      resid = resid-ctbasemake(sigma,bases,2)*c;
    end
    resid = reshape(resid,[n 1]);
  end
  v = reshape(v,[n 1]);
  x = reshape(x,[n size(x,2) 1]);
end