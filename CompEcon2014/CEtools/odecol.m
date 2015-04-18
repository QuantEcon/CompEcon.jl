%% ODECOL
%
%  Collocation ODE Boundary Value Problem Solver
%
%  ODECOL uses the method of collocation to solve a first-order ordinary
%  differential equation subject to boundary value constraints:
%    x'(t) = f(x(t))     t in [0,T]
%  s.t.
%    b_j(x(t_j)) = 0,   j=1..d
%  Here, x is d.1 vector-valued function defined on time domain [0,T] and
%  x' is its d.1 vector-valued derivative with respect to tome.
%
%  Usage
%    [x,t,resid,basis,c] = odecol(f,bv,T,n,bt,bx,type,c,varargin)
%  Let
%    d  = dimension of state process x
%    m  = number of time nodes at which approximation residual is returned
%  Input
%    f        : name of velocity function (see below)
%    bv       : d.1 boundary values 
%    T        ; time horizon
%    n        : number of collocation nodes on [0,T] (100)
%    bt       : d.1 boundary times (vector of zeros)
%    bx       : d.1 indices of variables in boundary conditions (1..d)
%    type     : basis function family for collocation - 'spli' or ('cheb')
%    c        : n.d guesses for basis function coefficients (array of zeros)
%    varargin : parameters passed to velocity function (empty)
%  Output
%    x        : m.d state process values at times t
%    t        : m.1 equally spaced time nodes on [0,T]
%    resid    : m.d residuals at times t
%    basis    : basis structured array
%    c        : n.d basis function coefficients
%    Note: The values of the state process and the residuals, if requested,
%    are returned on a grid of m equally-spaced time nodes t. The number m
%    of nodes in the grid is governed by nf, a user-supplied option with
%    default value of 10 that may be set by the user using optset (see
%    below).  If nf>0, m=nf*n; otherwise, the routine will return the state
%    process values and residuals on the n original time collocation nodes.
%  Velocity Function
%    User-supplied function that returns velocity at states according to
%    Format
%      v = f(x,varargin)
%    Input
%      x        : d.1 states
%      varargin : optional parameters passed to f
%    Output
%      v        : d.1 velocities
%  Options
%    nf        : grid expansion factor (10, see above)
%    showiters : display results of each iteration (0)

%  Copyright(c) 1997-2014
%    Mario J. Miranda - miranda.4@osu.edu
%    Paul L. Fackler  - paul_fackler@ncsu.edu
%
%  Uses BROYDEN

function [x,t,resid,basis,c] = odecol(f,bv,T,n,bt,bx,type,c,varargin)

fprintf('\nSolving ODE by collocation.\n')

% Set user options to defaults, if not set by user with OPTSET (see above)
nf        = optget('odecol','nf',10);
showiters = optget('odecol','showiters',1);
optset('broyden','showiters',showiters)

% Set inputs to default if nonexistent (see above)
if nargin<8, c    = []; end
if nargin<7, type = []; end
if nargin<6, bx   = []; end
if nargin<5, bt   = []; end
if nargin<4, n    = []; end
if nargin<3, error('Three Inputs Are Required'); end
if isempty(f),      error('Missing Velocity Function File'); end 
if isempty(bv),     error('Missing Boundary Values');        end  
if isempty(T),      error('Missing Time Horizon');           end  
if isempty(n),      n = 100;                                 end 
if isempty(bt),     bt = zeros(length(bv),1);                end 
if isempty(bx),     bx = (1:length(bv))';                    end 
if isempty(type),   type = 'cheb';                           end 
d = length(bv);

% Compute collocation nodes
if strcmp(type,'cheb')
  t = funnode(fundefn('cheb',n-1,0,T));
else
  t = funnode(fundefn('spli',n-1,0,T));
end

% Intialize coefficient vector if not provided
if isempty(c)
  if strcmp(type,'cheb')
    c = zeros(n,d);
    c(1,:) = bv';
  else
    c = ones(n,1)*bv';
  end
end

% Approximation Structure
basis = fundefn(type,n,0,T);

% Basis matrices
Phi  = funbase(basis,t);
Phi1 = funbase(basis,t,1);
phi  = funbase(basis,bt);

% Call rootfinding algorithm
optset('broyden','tol',1e-12)
c = broyden(@residual,c(:),Phi,Phi1,phi,f,bv,bx,varargin{:});
c = reshape(c,n,d);

% Compute solution at plotting nodes
if nf>0
  m = nf*n;
  t = nodeunif(m,0,T);
end
x = funeval(c,basis,t);

% Compute residual
if nargout>2
  dx = funeval(c,basis,t,1);
  resid = dx - feval(f,x',varargin{:})';
end


% Residual function used by odecol
function r = residual(c,Phi,Phi1,phi,f,bv,bx,varargin)

% Reshape coefficient vector
n = size(Phi,2);
d = length(bv);
c = reshape(c,n,d);

% Compute residuals at nodal times
x  = Phi*c;
dx = Phi1*c;
r  = dx - feval(f,x',varargin{:})';

% Compute residuals at boundaries
x = phi*c;
for j=1:d
  b = x(j,bx(j))-bv(j);
  r = [r(:);b];
end