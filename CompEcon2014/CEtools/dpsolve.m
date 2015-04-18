%% DPSOLVE
%
%  General Discrete-Time Bellman Equation Solver
%
%  Uses the method of collocation to solve finite- and infinite-horizon
%  discrete-time stochastic dynamic optimization models with discrete,
%  continuous, or mixed states and actions of arbitrary dimension.  Usage
%  of DPSOLVE is fully documented in the pdf file `DPSOLVE - A MATLAB
%  Routine for Solving Discrete-Time Bellman Equations' included in the
%  directory CompEcon2014 accompanying the distribution of this file.  A
%  summary of its features is provided here:
%
%  Usage
%    [c,sr,vr,xr,resid] = dpsolve(model,basis,v,x)
%  Let
%    ds = dimension of continuous state s
%    dx = dimension of continuous action x
%    de = dimension of continuous state transition shock e
%    ni = number of discrete states i
%    nj = number of discrete actions j
%    ne = number of continuous state transition shocks
%    n  = number of basis functions per continous state dimension (1.ds)
%    nc = number of nodes on continuous state collocation grid (prod(n))
%    ns = number of continuous state nodes output
%    nx = number of discretized continuous actions, if applicable
%  Input
%    model  : structured array containing model specifications (see below)
%    basis  : ds-dimensional basis for functions defined on continuous
%             state space
%    v      : nc.ni.nj initial guess for values at the continuous state
%             collocation nodes, per discrete state and discrete action
%             (optional, default is array of zeros)
%    x      : nc.dx.ni.nj initial guess for optimal continuous actions at
%             the continous state collocation nodes, per discrete state
%             and discrete action (optional, default is array of zeros)
%  Output
%    c      : nc.ni value function approximant basis function coefficients, 
%             per discrete state
%    sr     : ns.ds refined continuous state node grid
%    vr     : ns.ni.nj values at refined continuous state nodes, per 
%             discrete state and discrete action
%    xr     : ns.dx.ni.nj optimal continuous actions at refined continuous
%             state nodes, per discrete state and discrete action
%    resid  : ns.ni Bellman equation residuals at refined continuous state 
%             nodes, per discrete state
%    Note: If the residual is requested, then values, optimal actions, and
%    residuals are returned on a refined grid of ns continuous state nodes.
%    The degree of refinement is governed by nr, an optional parameter with
%    default value of 10 that may be set by the user using optset.  If
%    nr>0, the refined continuous state node grid is created by forming the
%    Cartesian product of nr*n equally-spaced coordinates along each
%    continous state dimension. If nr=0, the routine will return the
%    values, optimal actions, and residuals at the original collocation
%    continuous state nodes. On output, the singleton dimensions of c, vr,
%    xr, and resid are eliminated.
%  Model Structure
%    The structured array "model" contains different fields that specify
%    essential features of the model to be solved (default values for
%    optional fields in parentheses):
%      horizon   : time horizon (infinite)
%      func      : name of function file (required)
%      params    : model parameters required by function file (empty) 
%      discount  : discount factor (required)
%      ds        : dimension of continuous state (1)
%      dx        : dimension of continuous action (1)
%      ni        : number of discrete states (0)
%      nj        : number of discrete actions (0)
%      e         : ne.de discretized continuous state transition shocks (0)
%      w         : ne.1 continuous state shock probabilities (1)
%      q         : ni.ni.nj discrete state transition probabilities (empty)
%      h         : nj.ni deterministic discrete state transitions (empty)
%      X         : nx.dx discretized continuous actions (empty)
%    Note: If X is empty (the default), the routine will attempt to solve
%    the continuous action maximization problem embedded in the Bellman
%    equation by solving the associated Karush-Kuhn-Tucker conditions as a
%    nonlinear complementarity problem, using a derivative-based adapted
%    Newton method (see Miranda and Fackler).
%    Note: Either discrete state transition probabilities q or a
%    deterministic state transition function h, but not both, should be
%    specified by the user.
%  Function File
%    User-supplied routine that returns the bound, reward, and continuous
%    state transition functions and derivatives with respect to the
%    continuous action x at an arbitrary number ns continous states and
%    actions according to the format:
%      [out1,out2,out3] = func(flag,s,x,i,j,in,e,<params>)
%    Function File Input
%      flag      : flag that specifies function to be evaluated
%      s         : ns.ds continuous state nodes
%      x         : ns.dx continuous actions
%      i         : ns.1 or 1.1 discrete state indices between 1-ni
%      j         : ns.1 or 1.1 discrete action indices between 1-nj
%      in        : ns.1 or 1.1 next period discrete state indices between 1-ni
%      e         : ns.de continuous state transition shocks
%      params    : user-supplied list of function parameters
%    Function File Output
%      flag = 'b' returns lower and upper bounds on continuous action x
%        out1    : ns.dx lower bounds
%        out2    : ns.dx upper bounds
%        out3    : empty
%      flag = 'f' returns reward function values and derivatives
%        out1    : ns.1 reward function values
%        out2    : ns.dx first derivatives with respect to x
%        out3    : ns.dx.dx second derivatives with respect to x
%      flag = 'g' returns transition function values and derivatives
%        out1    : ns.ds transition function values
%        out2    : ns.ds.dx first derivatives with respect to x
%        out3    : ns.ds.dx.dx second derivatives with respect to x
%      Note: If the continuous action is discretized, then the function
%      file need only provide the reward and transition values, with the
%      reward set to negative infinity at non-feasible values of x and j.
%  Options
%    algorithm : collocation equation solution algorithm, either Newton's 
%                method (`newton') or function iteration `funcit'
%    ncpmethod : nonlinear complementarity solution algorithm for 
%                continuous action maximization problem embedded in Bellman
%                equation, either semi-smooth formulation 'ssmooth' or
%                min-max formulation 'minmax', the default
%    maxit     : maximum number of iterations, collocation equation 
%                solution algorithm (500)
%    maxitncp  : maximum number of iterations, nonlinear complementarity 
%                solver (50)
%    tol       : convergence tolerance (square root of machine epsilon)
%    nr        : continuous state grid refinement factor (10)
%    output    : whether to generate printed output (1)

%  Copyright(c) 1997-2014
%    Mario J. Miranda - miranda.4@osu.edu
%    Paul L. Fackler  - paul_fackler@ncsu.edu

function [c,sr,vr,xr,resid] = dpsolve(model,basis,v,x)

% Set user options to defaults, if not set by user with OPTSET (see above)
algorithm = optget('dpsolve','algorithm','newton');
tol       = optget('dpsolve','tol',sqrt(eps));
maxit     = optget('dpsolve','maxit',500);  
nr        = optget('dpsolve','nr',10);  
output    = optget('dpsolve','output',1);  
if nargout<5
  nr = 0;
end

% Set model fields to default if nonexistent (see above)
if ~isfield(model,'horizon'),  model.horizon = inf;              end  
if ~isfield(model,'func'),     error('Missing Function File');   end  
if ~isfield(model,'params'),   error('Missing Parameter List');  end  
if ~isfield(model,'discount'), error('Missing Discount Factor'); end  
if ~isfield(model,'ds'),       model.ds = 1;                     end  
if ~isfield(model,'dx'),       model.dx = 1;                     end  
if ~isfield(model,'ni'),       model.ni = 1;                     end  
if ~isfield(model,'nj'),       model.nj = 1;                     end  
if ~isfield(model,'e'),        model.e  = 0;                     end  
if ~isfield(model,'w'),        model.w  = 1;                     end  
if ~isfield(model,'q'),        model.q  = [];                    end
if ~isfield(model,'h'),        model.h  = [];                    end
if ~isfield(model,'X'),        model.X  = zeros(1,model.dx);     end
model.ni = max(1,model.ni);
model.nj = max(1,model.nj);

% Unpack model structure
T      = model.horizon;  
func   = model.func;  
params = model.params;
ds     = model.ds;  
dx     = model.dx;  
ni     = model.ni;  
nj     = model.nj;  
e      = model.e;  
w      = model.w;  
q      = model.q;
h      = model.h;

% Equivalent transition probability matrix for deterministic discrete state
if ni==1; h = ones(nj,1); end
if isempty(q)
  q = zeros(ni,ni,nj);
  if isempty(h)
    if ni==nj
      disp('Neither q nor h specified; will default to h(j,i)=j.')
      disp('Please specify q or h if this is not correct.')
      disp(' ')
      h = (1:nj)'*ones(1,ni);
    else
      error('Either q or h must be specified.')
    end
  end
  for i=1:ni
    for j=1:nj
      q(i,h(j,i),j) = 1;
    end
  end
end
model.q = q;
model.h = [];

% Determine number of continuous state collocation nodes
n  = basis.n;               % number of continuous state collocation nodes
                            % ... per dimension (1.ds)
nc = prod(n);               % number of continuous state collocation nodes
  
% Compute continuous state collocation nodes and collocation matrix
s   = funnode(basis);      	% continuous state collocation nodes
Phi = funbase(basis,s);     % collocation matrix

% Initial guess for values and actions if not provided
if nargin<3
  v = zeros(nc,ni,nj);
else
  if T<inf
    v = reshape(v,nc,ni,1);
    v = v(:,:,ones(1,nj));
  else
    v = reshape(v,nc,ni,nj);
  end
end
if nargin<4, x = zeros(nc,dx,ni,nj); else x = reshape(x,nc,dx,ni,nj); end

% Initial or post-terminal value function approximant basis coefficients
c = Phi\max(v,[],3);

% Solve collocation equations for finite horizon model by backward recursion
if T<inf
  if output, display('Solving finite-horizon model by backward recursion.'), end
  cc = zeros(nc,ni,T+1);
  vv = zeros(nc,ni,nj,T+1);
  xx = zeros(nc,dx,ni,nj,T);
  cc(:,:,T+1) = c;
  vv(:,:,:,T+1) = v;
  for t=T:-1:1
    if output, fprintf ('Period %5i\n',t), end  % print period
    [v,x] = vmax(model,basis,s,x,c);            % update optimal values and actions
    c = Phi\max(v,[],3);                        % update basis coefficients
    cc(:,:,t) = c;
    vv(:,:,:,t) = v;
    xx(:,:,:,:,t) = x;
  end
  sr = s;
  c  = squeeze(cc);
  vr = squeeze(vv);
  xr = squeeze(xx);
  resid = [];
  return
end

tic

% Solve collocation equation for infinite-horizon model by function iteration.
if T==inf && strcmp(algorithm,'funcit')
  if output, display('Solving infinite-horizon model collocation equation by function iteration.'), end
  for it=1:maxit
    cold = c;                                           % store old basis coefficients
    [v,x] = vmax(model,basis,s,x,c);                    % update optimal values and actions
    c = Phi\max(v,[],3);                                % update basis coefficients
    change = norm(c(:)-cold(:),inf);                    % compute change
    if output, fprintf ('%4i %10.1e\n',it,change), end	% print change
    if change<tol, break, end;                          % convergence test
  end
end

% Solve collocation equation for infinite horizon model by Newton's method.
if T==inf && strcmp(algorithm,'newton')
  if output, display('Solving infinite horizon model collocation equation by Newton''s method.'), end
  Phik = kron(eye(ni),Phi);
  for it=1:maxit
    cold = c;                                           % store old basis coefficients
    [v,x,vc] = vmax(model,basis,s,x,c);                 % update optimal contingent values and actions
    vm = max(v,[],3);                                   % compute optimal value
    c = c(:)-(Phik-vc)\(Phik*c(:)-vm(:));               % update basis coefficients
    c = reshape(c,nc,ni);                               % reshape basis coefficient array
    change = norm(c(:)-cold(:),inf);                    % compute change
    if output, fprintf ('%4i %10.1e\n',it,change), end	% print change
    if change<tol, break, end;                          % convergence test
  end
end

if output
  if it==maxit, display('Failure to converge in dpsolve.'), end
  fprintf('Elapsed Time = %7.2f Seconds\n',toc)
end

% Check whether continuous state transitions remain within approximation interval
snmin =  inf;
snmax = -inf;
for i=1:ni
  [~,jmax] = max(v(:,i,:),[],3);
  for in=1:ni
    if q(i,in)==0, continue, end
    for j=1:nj
      is = find(jmax==j);
      if isempty(is), continue, end
      ns = length(is);
      for k=1:size(e,1);
        ee = e(k+zeros(ns,1),:);
        sn = feval(func,'g',s(is,:),x(is,:,i,j),i,j,in,ee,params{:});
        snmin = min(snmin,min(sn));
        snmax = max(snmax,max(sn));
      end
    end
  end
end
if output
  if any(snmin<basis.a-eps), display('Extrapolating below smin. Decrease smin.'), end;
  if any(snmax>basis.b+eps), display('Extrapolating above smax. Increase smax.'), end;
  disp([basis.a' snmin' snmax' basis.b'])
end

% Compute values and actions on refined continuous state grid, if requested.
if nr>0
  n  = nr*n;
  ns = prod(n);
  for is=1:ds
    srcoord{is} = linspace(basis.a(is),basis.b(is),n(is))';
  end
  sr = gridmake(srcoord);
  if dx>0
    xr = zeros(ns,dx,ni,nj);
    for i=1:ni
      for j=1:nj
        xr(:,:,i,j) = funeval(Phi\x(:,:,i,j),basis,sr);
      end
    end
  else
    xr = [];
  end
  [vr,xr] = vmax(model,basis,sr,xr,c);
else
  sr = s;
  vr = v;
  xr = x;
end

% Compute Bellman equation residual on refined continuous state grid, if requested
if nargout>4
  vrproj = funeval(c,basis,sr);
  resid = vrproj - max(vr,[],3);
  resid = squeeze(resid);
end

% Eliminate singleton dimensions to facilitate analysis in calling program
c  = squeeze(c);
vr = squeeze(vr);
xr = squeeze(xr);

  
%% VMAX
%
%  Routine called by dpsolve to maximize the optimand embedded in the
%  Bellman equation with respect to the continuous action x at ns input
%  continuous state nodes s and discrete states i, coningent on specified
%  discrete actions j.  Does so by solving associated KKT conditions as a
%  nonlinear complementarity problem, using a derivative-based adapted
%  Newton method (see Miranda and Fackler). Routine optionally also
%  computes derivative of the optimand with respect to the value function
%  approximant basis function coefficients, which is required to solve the
%  collocation equation using Newton's method.
%
%  Usage
%   [v,x,vc] = vmax(model,basis,s,x,c)
%  Let
%    ds = dimension of continuous state
%    dx = dimension of continuous action
%    de = dimension of continuous state transition shock
%    ni = number of discrete states
%    nj = number of discrete actions
%    ne = number of continuous state transition shocks
%    n  = number of basis functions per continous state dimension (1.ds)
%    nc = number of nodes on continuous state collocation grid (prod(n))
%    ns = number of continuous state nodes input
%    nx = number of discretized continuous actions, if applicable
% Input
%   model : structured array containing model specifications 
%   basis : ds-dimensional basis for functions defined on continuous
%            state space
%   s     : ns.ds continuous state node grid
%   x     : nx.dx.ni.nj initial guess for optimal continuous actions at
%           the continous state collocation nodes, per discrete state
%           and discrete action
%    c    : nc.ni value function approximant basis function coefficients, 
%           per discrete state
% Output
%   v     : ns.1 values at state nodes
%   x     : ns.dx optimal continuous actions at state nodes
%   vc    : nn.nn (nn=nc*ni*nj) array of derivatives of values with
%           respect to basis coefficients

function [v,x,vc] = vmax(model,basis,s,x,c)

% Unpack user options
tol       = optget('dpsolve','tol',sqrt(eps));
maxitncp  = optget('dpsolve','maxitncp',50);  
ncpmethod = optget('dpsolve','ncpmethod','minmax'); 

% Unpack model structure
func   = model.func;  
params = model.params;
delta  = model.discount;
ds     = model.ds;  
dx     = model.dx;  
ni     = model.ni;  
nj     = model.nj;  
e      = model.e;  
w      = model.w;  
q      = model.q;
X      = model.X;

% Determine number of continuous state nodes
ns = size(s,1);             % number of continuous state nodes
nx = size(X,1);             % number of discretized continuous actions

% Compute Bellman equation optimand - discrete or discretized continuous action
if dx==0||nx>1
  v = zeros(ns,ni,nj);
  x = zeros(ns,dx,ni,nj);
  for i=1:ni
    for j=1:nj
      vv = zeros(ns,nx);
      if nx>1
        [xl,xu] = feval(func,'b',s,[],i,j,[],[],params{:});
      end
      for ix=1:nx
        vv(:,ix) = -inf;
        xx = X(ix+zeros(ns,1),:);
        if nx>1
          is = find(all(xx>=xl,2).* all(xx<=xu,2)==1);
        else
          is = 1:ns;
        end
        if ~isempty(is)
          % Initialize optimand and derivatives with reward terms
          vv(is,ix) = feval(func,'f',s(is,:),xx(is,:),i,j,[],[],params{:});
          for k=1:length(w)
            % Compute states next period and derivatives
            ee = e(k+zeros(length(is),1),:);
            for in=1:ni
              if q(i,in,j)==0, continue, end
              snext = feval(func,'g',s(is,:),xx(is,:),i,j,in,ee,params{:});
              snext = real(snext);
              prob = w(k)*q(i,in,j);
              vn = funeval(c(:,in),basis,snext);
              vv(is,ix) = vv(is,ix) + delta*prob*vn;
            end
          end
        end
      end
      [v(:,i,j),ix] = max(vv,[],2);
      x(:,:,i,j) = X(ix,:);
    end
  end
end

% Compute Bellman equation optimand - continuous or mixed action
if dx>0&&nx<2
  x = reshape(x,ns,dx,ni,nj);
  v = zeros(ns,ni,nj);
  for i=1:ni
    for j=1:nj
      [xl,xu] = feval(func,'b',s,[],i,j,[],[],params{:});
      for it=1:maxitncp
        % Initialize optimand and derivatives with reward terms
        [vv,vx,vxx] = feval(func,'f',s,x(:,:,i,j),i,j,[],[],params{:});
        for k=1:length(w)
          % Compute states next period and derivatives
          ee = e(k+zeros(ns,1),:);
          for in=1:ni
            if q(i,in,j)==0, continue, end
            [snext,snx,snxx] = feval(func,'g',s,x(:,:,i,j),i,j,in,ee,params{:});
            snext = real(snext);
            B = funbasex(basis,snext,[0;1;2]);
            vnss = zeros(ns,nj,ds,ds);
            prob = w(k)*q(i,in,j);
            vn  = funeval(c(:,in),basis,B);
            vns = funeval(c(:,in),basis,B,eye(ds));
            for is=1:ds
              for js=is:ds
                order = zeros(1,ds);
                order(is) = order(is)+1;
                order(js) = order(js)+1;
                vnss(:,is,js) = funeval(c(:,in),basis,B,order);
                vnss(:,js,is) = vnss(:,is,js);
              end
            end
            vv = vv + delta*prob*vn;
            for ix=1:dx
              for is=1:ds
                vx(:,ix) = vx(:,ix) + delta*prob*vns(:,is).*snx(:,is,ix);
                for jx=1:dx
                  vxx(:,ix,jx) = vxx(:,ix,jx) + delta*prob*vns(:,is).*snxx(:,is,ix,jx);
                  for js=1:ds
                    vxx(:,ix,jx) = vxx(:,ix,jx) + delta*prob*vnss(:,is,js).*snx(:,is,ix).*snx(:,js,jx);
                  end
                end
              end
            end
          end
        end
        % Compute Newton step, update continuous action, check convergence
        [vx,delx] = lcpstep(ncpmethod,x(:,:,i,j),xl,xu,vx,vxx);
        x(:,:,i,j) = x(:,:,i,j) + delx;
        if norm(vx(:),inf)<tol, break, end;
      end
      v(:,i,j) = vv;
    end
  end
end

if nargout<3, return, end

% Compute derivative of Bellman equation optimand with respect to basis
% coefficients for Newton method, if requested

if ni*nj>1
  vc = zeros(ns,ni,ns,ni);
  [~,jmax] = max(v,[],3);
  for i=1:ni
    for j=1:nj
      is = find(jmax(:,i)==j);
      if isempty(is), continue, end
      for k=1:length(w)
        ee = e(k+zeros(ns,1),:);
        for in=1:ni
          if q(i,in,j)==0, continue, end
          snext = feval(func,'g',s(is,:),x(is,:,i,j),i,j,in,ee(is,:),params{:});
          B = full(funbase(basis,snext));
          prob = w(k)*q(i,in,j);
          vc(is,i,:,in) = vc(is,i,:,in) + delta*prob*reshape(B,length(is),1,ns);
        end
      end
    end
  end
  vc = reshape(vc,ni*ns,ni*ns);
else
  vc = zeros(ns,ns);
  for k=1:length(w)
    ee = e(k+zeros(ns,1),:);
    snext = feval(func,'g',s,x,1,1,1,ee,params{:});
    vc = vc + delta*w(k)*funbase(basis,snext);
  end
end