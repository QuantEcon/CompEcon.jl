%% GAMESOLVE
%
%  Discrete-Time Continuous-State Continuous-Action Infinite-Horizon 
%  Multi-Agent Dynamic Game Model Bellman Equation Solver
%
%  Uses the method of collocation to compute an approximate solution to
%  the Bellman equation of a discrete-time infinite-horizon multi-agent
%  dynamic game model with continuous state variable s and continuous 
%  joint action variable x. Returns value function basis function 
%  coefficients c and optimal values v, actions x, and Bellman equation 
%  residuals resid on a refined grid of state nodes.  Specifically, solves
%     V_i(s) = max_x_i {f_i(s,x_i,x_(i))+delta*E_eV(g(s,x_i,x_(i),e))}
%  for i=1..dp.
%
%  Usage
%    [c,sr,vr,xr,resid] = gamesolve(model,basis,v,x)
%  Let
%    nc = number of state collocation nodes
%    ns = number of state nodes on refined output grid
%    dp = number of players
%    ds = dimension of state s
%    dx = dimension of individual action x_i
%    de = dimension of shock e
%  Input
%    model : model structure
%    basis : approximation structure
%    v     : nc.dp     initial guess for optimal values at collocation nodes
%    x     : nc.dx.dp  initial guess for optimal actions at collocation nodes
%  Output
%    c     : nc.dp     value function approximant basis function coefficients
%    sr    : ns.ds     refined state nodes
%    vr    : ns.dp     optimal values at refined state nodes
%    xr    : ns.dx.dp  optimal actions at refined state nodes
%    resid : ns.dp     Bellman equation residuals at refined state nodes
%    Note: The number of refined state nodes ns is determined by the grid 
%    refinement factor nres, an optional parameter with default value of 10
%    that may be set by the user. If nres>0, the routine will return values 
%    and residuals on a refined grid of ns uniformly-spaced state nodes that 
%    possesses, along each dimension, nres times the coordinates possessed
%    by the state collocation node grid. If nres=0 is set by the user, the 
%    routine will return the values and residuals on the original state 
%    collocation nodes.  
% 
%  Model Structure Fields
%    func      : name of function file (see below)
%    params    : parameters passed to function file
%    discount  : discount factor
%    dp        : number of players
%    ds        : dimension of state s
%    dx        : dimension of individual action x_i
%    e         : ne.de shocks
%    w         : ne.1  shock probabilities 
% 
%  Function File
%    A user-supplied function that returns the bound, reward, and state 
%    transitions and their first and second derivatives with respect to 
%    player i's action x_i, at an arbitrary number ns of states s and joint 
%    actions x, according to the format
%      [out1,out2,out3] = func(flag,i,s,x,e,params)
%    Function File Input
%      flag      : flag indicating function to be evaluated 
%      i         : player index, an integer between 1 and dp
%      s         : ns.ds    states
%      x         : ns.dx.dp actions
%      e         : ns.de    shocks
%      params    : parameters passed to function file
%    Function File Output
%      if flag = 'b', returns bounds on individual action x_i
%        out1    : ns.dx lower bounds on individual action x_i
%        out2    : ns.dx upper bounds on individual action x_i
%        out3    : empty
%      if flag = 'f', returns reward and derivatives w.r.t. x_i
%        out1    : ns.1     reward f_i
%        out2    : ns.dx    first derivative of f_i with respect to x_i
%        out3    : ns.dx.dx second derivative of f_i with respect to x_i
%      if flag = 'g', returns state transition and derivatives w.r.t. x_i
%        out1    : ns.ds       state g
%        out2    : ns.ds.dx    first derivative of g with respect to x_i
%        out3    : ns.ds.dx.dx second derivative of g with respect to x_i
% 
%  Options
%    These options may be set by user with OPTSET (defaults in parentheses):
%    algorithm : collocation equation solution algorithm, either Newton
%                method ('newton') or function iteration 'funcit'
%    tol       : solution algorithm convergence tolerance (sqrt(eps))
%    maxit     : solution algorithm maximum number of iterations (500)
%    nres      : state grid refinement factor (10)

% Copyright(c) 1997-2010
%    Mario J. Miranda - miranda.4@osu.edu
%    Paul L. Fackler  - paul_fackler@ncsu.edu

function [c,sr,vr,xr,resid] = gamesolve(model,basis,v,x)

% Set options to defaults, if not set by user with OPTSET (see above)
algorithm = optget('gamesolve','algorithm','newton');
tol       = optget('gamesolve','tol',sqrt(eps));
maxit     = optget('gamesolve','maxit',500);
nres      = optget('gamesolve','nres',10);
limcheck  = optget('gamesolve','limcheck',1);

% Set model structure parameters to defaults, if nonexistent
if ~isfield(model,'dp'), model.dp = 2; end % number of players
if ~isfield(model,'ds'), model.ds = 1; end % dimension of state s
if ~isfield(model,'dx'), model.dx = 1; end % dimension of individual action x_i
if ~isfield(model,'e'),  model.e  = 0; end % shocks
if ~isfield(model,'w'),  model.w  = 1; end % shock probabilities

% Unpack model structure
func   = model.func;        % name of function file
params = model.params;      % parameters passed to function file
dp = model.dp;              % number of players
ds = model.ds;              % dimension of state s
dx = model.dx;              % dimension of individual action x_i
e  = model.e;               % shocks

% Determine number of collocation nodes
n  = basis.n;               % number of state collocation node coordinates by dimension
nc = prod(n);               % number of state collocation nodes

% Set initial guess for values and actions if not provided
if nargin<3, v = zeros(nc,dp); end
if nargin<4, x = zeros(nc,dx,dp); end
x = reshape(x,nc,dx,dp);
  
% Compute collocation nodes and matrix and basis coefficients
s = funnode(basis);         % collocaton nodes
Phi = funbase(basis,s);     % collocation matrix
c = Phi\v; 

% Solve collocation equation
tic
switch algorithm
  case 'funcit'
    display('Solving collocation equation by function iteration.')
    for it=1:maxit
      cold = c;                             	% store old basis coefficients
      [v,x] = vmax(model,basis,s,x,c);        % update optimal values and actions
      c = Phi\v;                              % update basis coefficient
      change = norm(c(:)-cold(:),inf);       	% compute change
      fprintf ('%4i %10.1e\n',it,change)      % print change
      if change<tol, break, end;              % convergence check
    end
  case 'newton'
    display('Solving collocation equation by Newton method.')
    for it=1:maxit
      cold = c;                               % store old basis coefficients
      [v,x,vc] = vmax(model,basis,s,x,c);     % update optimal values and actions
      for ip=1:dp                             % update basis coefficients
        c(:,ip) = c(:,ip)-(Phi-vc(:,:,ip))\(Phi*c(:,ip)-v(:,ip));    
      end
      change = norm(c(:)-cold(:),inf);      	% compute change
      fprintf ('%4i %10.1e\n',it,change)      % print change
      if change<tol, break, end;              % convergence check
    end
end
if it==maxit, display('Failure to converge in gamesolve.'), end
fprintf('Elapsed Time = %7.2f Seconds\n',toc)

% Check whether state transitions remain within approximation interval
if limcheck
  snmin =  inf;
  snmax = -inf;
  for k=1:size(e,1);
    ee = e(k*ones(nc,1),:);
    sn = feval(func,'g',1,s,x,ee,params{:});
    snmin = min(snmin,min(sn));
    snmax = max(snmax,max(sn));
  end
  if any(snmin<basis.a-eps), display('Extrapolating below smin. Decrease smin.'), end;
  if any(snmax>basis.b+eps), display('Extrapolating above smax. Increase smax.'), end;
  disp([basis.a' snmin' snmax' basis.b'])
end
  
% Compute values and actions on refined state grid, if requested.
if nres>0
  nr = nres*n;
  for is=1:ds
    srcoord{is} = linspace(basis.a(is),basis.b(is),nr(is))';
  end
  sr = gridmake(srcoord);
  ns = prod(nr);
  xr = zeros(ns,dx,dp);
  for i=1:dp
    xr(:,:,i) = funeval(Phi\x(:,:,i),basis,sr);
  end
  [vr,xr] = vmax(model,basis,sr,xr,c);
else
  sr = s;
  vr = v;
  xr = x;
end

% Compute Bellman equation residual on output state grid
resid = funeval(c,basis,sr) - vr;


%% VMAX
%
%   Maximizes optimand of a discrete-time infinite-horizon Bellman equation 
%   with respect to the continuous action x at ns continuous state nodes s
%   by solving the Karush-Kuhn-Tucker necessary conditions as a nonlinear 
%   complementarity problem.  Uses successive linearization, which requires, 
%   at each iteration, evaluation of Bellman optimand and its derivatives 
%   with respect to action x:
%       v   = f   + delta*E(vn(g))
%       vx  = fx  + delta*E(vns*gx)
%       vxx = fxx + delta*E(vns*gxx+vnss*gx*gx)
%   Optionally, will also compute derivative of the optimand with respect
%   to the value function approximant basis function coefficients. 
%
% Usage
%   [v,x,vc] = vmax(model,basis,s,x,c)     Called by gamesolve
% Let
%   nc = number of state collocation nodes
%   ns = number of state nodes input
%   dp = number of players
%   ds = dimension of state s
%   dx = dimension of individual action x_i
% Input
%   model : model structure
%   basis : approximation structure
%   s     : ns.ds     state nodes
%   x     : ns.dx.dp  guess for optimal actions at state nodes
%   c     : nc.dp     value function approximant basis function coefficients
% Output
%   v     : ns.dp     optimal values at state nodes
%   x     : ns.dx.dp  optimal actions at state nodes
%   vc    : nc.nc     derivatives of values with respect to basis coefficients
% Options
%   These options may be set by user with OPTSET (defaults in parentheses):
%   tol       : NCP algorithm convergence tolerance (sqrt(eps))
%   maxit     : NCP algorithm maximum number of iterations (50)
%   lcpmethod : linear complementarity solution algorithm, either mini-max
%               function ('minmax') or Fischer function 'ssmooth' 

function [v,x,vc] = vmax(model,basis,s,x,c)

% Set options to defaults, if not set by user with OPTSET (see above)
tol       = optget('vmax','tol',sqrt(eps));
maxit     = optget('vmax','maxit',50);
lcpmethod = optget('vmax','lcpmethod','minmax');

% Unpack model structure
func   = model.func;        % name of function file
params = model.params;      % parameters passed to function file
delta  = model.discount;    % discount factor
dp     = model.dp;          % number of players
ds     = model.ds;          % dimension of state s
dx     = model.dx;          % dimension of action x
e      = model.e;           % shocks
w      = model.w;           % shock probabilties

% Determine number of state nodes
ns = size(s,1);

% Intialize arrays
v = zeros(ns,dp);

% Compute Bellman equation optimand
for i=1:dp
  [xl,xu] = feval(func,'b',i,s,x,[],params{:});
  for it=1:maxit
    % Initialize optimand and derivatives
    [vi,vx,vxx] = feval(func,'f',i,s,x,[],params{:});
    vx  = reshape(vx,ns,1,dx);
    vxx = reshape(vxx,ns,dx,dx);
    for k=1:length(w)
      % Compute states next period and derivatives
      ee = e(k*ones(ns,1),:);
      [snext,snx,snxx] = feval(func,'g',i,s,x,ee,params{:});
      B = funbasex(basis,snext,[0;1;2]);
      vn   = funeval(c(:,i),basis,B);
      vns  = funeval(c(:,i),basis,B,eye(ds));
      vnss = zeros(ns,ds,ds);
      for is=1:ds
        for js=is:ds
          order = zeros(1,ds);
          order(is) = order(is)+1;
          order(js) = order(js)+1;
          vnss(:,is,js) = funeval(c(:,i),basis,B,order);
          vnss(:,js,is) = vnss(:,is,js);
        end
      end
      vi = vi + delta*w(k)*vn;
      for ix=1:dx
        for is=1:ds
          vx(:,ix) = vx(:,ix) + delta*w(k)*vns(:,is).*snx(:,is,ix);
          for jx=1:dx
            vxx(:,ix,jx) = vxx(:,ix,jx) + delta*w(k)*vns(:,is).*snxx(:,is,ix,jx);
            for js=1:ds
              vxx(:,ix,jx) = vxx(:,ix,jx) + delta*w(k)*vnss(:,is,js).*snx(:,is,ix).*snx(:,js,jx);
            end
          end
        end
      end
    end
    % Compute Newton step, update action, check convergence
    v(:,i) = vi;
    [vx,deltax] = lcpstep(lcpmethod,x(:,:,i),xl,xu,vx,vxx);
    err = max(abs(vx),[],2);
    if all(err<tol), break, end;
    x(:,:,i) = x(:,:,i)+deltax;
  end
end

% Compute derivative of Bellman equation optimand with respect to basis
% coefficients for Newton method, if requested
if nargout>2
  vc = zeros(ns,ns,dp);
  for ip=1:dp
    for k=1:length(model.w)
      ee = e(k*ones(ns,1),:);
      g = feval(func,'g',ip,s,x,ee,params{:});
      phinext = funbase(basis,g);
      vc(:,:,ip) = vc(:,:,ip) + delta*w(k)*phinext;
    end
  end
end