%% DPSOLVE2
%
%  General Discrete-Time Bellman Equation Solver
%
%  DPSOLVE2 uses the method of collocation to solve finite- and infinite-
%  horizon discrete-time stochastic dynamic optimization models with
%  discrete, continuous, or mixed states and actions of arbitrary
%  dimensions.  DPSOLVE2 directly computes approximations for the
%  discrete-action-contingent value functions, whereas DPSOLVE computes
%  approximation for the value function.  If the model admits multiple
%  discrete actions, dpsolve2 is slower and requires more memory than
%  DPSOLVE, but can sometimes produce more accurate results. If the model
%  does not admit discrete actions, the two routines produce identical
%  results and do not differ with respect to memory usage or execution
%  time.  Usage of DPSOLVE2 is fully documented in the pdf file `DPSOLVE -
%  A MATLAB Routine for Solving Discrete-Time Bellman Equations' included
%  in the directory CompEcon2010 accompanying the distribution of this
%  file.  A summary of its features is provided here:
%
%  Usage
%    [c,sr,vr,xr,resid] = dpsolve2(model,basis,v,x)
%  Let
%    n  = number of basis functions and continuous state collocation nodes
%    ns = number of continuous state nodes on refined output array
%    ds = dimension of continuous state s
%    dx = dimension of continuous action x
%    de = dimension of continuous state transition shock e
%    ni = number of discrete states i
%    nj = number of discrete actions j
%    nx = number of discretized continuous actions, if applicable
%  Input
%    model  : structured array containing model specifications (see below)
%    basis  : n-dimensional basis for real-valued functions defined on the 
%             continuous state space
%    v      : n.ni.nj array of initial guesses for discrete-action-
%             contingent value function values at the n continuous state
%             collocation nodes, per discrete state and discrete action
%             (optional, default is array of zeros)
%    x      : n.dx.ni.nj array of initial guesses for optimal continuous
%             actions at the n state collocation nodes, per discrete state
%             and discrete action (optional, default is array of zeros)
%  Output
%    c      : n.ni.nj array of discrete-action contingent value function 
%             approximant basis function coefficients, per discrete state
%             and discrete action 
%    sr     : ns.ds refined array of continuous state nodes
%    vr     : ns.ni.nj array of discrete-action-contingent values on the
%             refined array of continuous state nodes, per discrete state
%             and discrete action
%    xr     : ns.dx.ni.nj array of discrete-action-contingent optimal
%             continuous actions on the refined array of continuous state
%             nodes, per discrete state and discrete action
%    resid  : ns.ni array of Bellman equation residuals on the refined 
%             array of continuous state nodes, per discrete state
%    Note: If the residual is requested, values, optimal continuous
%    actions, and residuals are returned on a refined array of ns
%    equally-spaced continuous state nodes. The degree of refinement is
%    governed by nr, an optional parameter with default value of 10 that
%    may be set by the user using optset (see below).  If nr>0, the refined
%    continuous state array is created by forming the Cartesian product of
%    equally-spaced coordinates along each dimension, with nr times the
%    number of coordinates possessed by the continuous state collocation
%    array along each dimension. If nr=0, the routine will return the
%    values, optimal continuous actions, and residuals on the original
%    array of continuous state collocation nodes.
%    Note: The singleton dimensions of c, vr, xr, and resid are eliminated
%    to facilitate analysis in the calling program.
%  Model Structure
%    The structured array "model" contains different fields that specify
%    essential features of the model to be solved (default values for
%    optional fields in parentheses):
%      horizon   : time horizon (infinite)
%      func      : name of function file (required)
%      params    : model parameters required by function file (empty) 
%      discount  : discount factor (required)
%      ds        : dimension ds of the continuous state s (1)
%      dx        : dimension dx of the continuous action x (0)
%      ni        : number ni of discrete states i (none)
%      nj        : number nj of discrete actions j (none)
%      e         : ne.de array of discretized continuous state shocks (0)
%      w         : ne.1 vector of discretized continuous state shock
%                  probabilities (1)
%      q         : ni.ni.nj array of stochastic discrete state transition
%                  probabilities (empty)
%      h         : nj.ni array of deterministic discrete state transitions
%                  (empty)
%      X         : nx.dx array of discretized continuous actions (empty)
%    Note: If X is empty (the default), the routine will attempt to solve
%    the continuous action maximization problem embedded in the Bellman
%    equation by solving the associated Karush-Kuhn-Tucker conditions as a
%    nonlinear complementarity problem, using a derivative-based adapted
%    Newton method (see Miranda and Fackler).
%    Note: If the stochastic discrete state transition probabilities q are
%    specified, then the deterministic state transitions h should not be
%    specified, and vice versa. 
%  Function File
%    User-supplied function that returns the bound, reward, and continuous
%    state transition functions and needed derivatives with respect to the
%    continuous action x at an arbitrary number ns of continuous states s
%    and continuous actions x, given a discrete state i and a discrete
%    action j, according to the format
%      [out1,out2,out3] = func(flag,s,x,i,j,e,<params>)
%    Function File Input
%      flag      : flag that specifies the desired function
%      s         : ns.ds array of continuous state nodes
%      x         : ns.dx array of continuous actions
%      i         : a single discrete state index, an integer between 1-ni
%      j         : a single discrete action index, an integer between 1-nj
%      e         : ns.de array of continuous state transition shocks
%      params    : user-supplied list of function parameters
%    Function File Output
%      flag = 'b' returns lower and upper bounds on continuous action x
%        out1    : ns.dx lower bounds on continuous action x
%        out2    : ns.dx upper bounds on continuous action x
%        out3    : empty
%      flag = 'f' returns reward and derivatives
%        out1    : ns.1 reward function f values
%        out2    : ns.dx first derivative of f with respect to x
%        out3    : ns.dx.dx second derivative of f with respect to x
%      flag = 'g' returns continuous state transition and derivatives
%        out1    : ns.ds continuous state transition function g values
%        out2    : ns.ds.dx first derivative of g with respect to x
%        out3    : ns.ds.dx.dx second derivative of g with respect to x
%      Note: If the continuous action is discretized, then the function
%      file need only provide the reward and transition, with the reward
%      set to negative infinity at non-feasible values of x as in
%      flag = 'f' returns reward and derivatives
%        out     : ns.1 reward f
%      flag = 'g' returns continuous state transition and derivatives
%        out     : ns.ds continuous state transition g
%  Options
%    algorithm : collocation equation solution algorithm, either Newton's 
%                method (`newton') or function iteration `funcit'
%    ncpmethod : nonlinear complementarity solution algorithm for 
%                continuous action maximization problem embedded in Bellman
%                equation, either semi-smooth formulation ('ssmooth') or
%                min-max formulation 'minmax'
%    maxit     : maximum number of iterations, collocation equation 
%                solution algorithm (500)
%    maxitncp  : maximum number of iterations, nonlinear complementarity 
%                solver (50)
%    tol       : convergence tolerance (square root of machine epsilon)
%    nr        : continuous state array refinement factor (10)
%    output    : whether to generate printed output (1)

%  Copyright(c) 1997-2011
%    Mario J. Miranda - miranda.4@osu.edu
%    Paul L. Fackler  - paul_fackler@ncsu.edu

function [c,sr,vr,xr,resid,ierr] = dpsolve2(model,basis,v,x)

% Set user options to defaults, if not set by user with OPTSET (see above)
algorithm = optget('dpsolve2','algorithm','newton');
tol       = optget('dpsolve2','tol',sqrt(eps));
maxit     = optget('dpsolve2','maxit',500);  
nr        = optget('dpsolve2','nr',10);  
output    = optget('dpsolve2','output',1);  
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
q      = model.q;
h      = model.h;

% Equivalent transition probability matrix for deterministic discrete state
if ni==1; h = ones(nj,1); end
if isempty(q)
  if isempty(h)
    if ni==nj
      'Neither q nor h specified; will default to h(j,i)=j.'
      'Please specify q or h if this is not correct.'
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
n  = basis.n;               % number of continuous state collocation
                            % node coordinates by dimension
nc = prod(n);               % number of continuous state collocation nodes

% Set initial guess for values and actions if not provided
if nargin<3, v = zeros(nc,ni,nj);    else v = reshape(v,nc,ni,nj); end
if nargin<4, x = zeros(nc,dx,ni,nj); else x = reshape(x,nc,dx,ni,nj); end
  
% Compute continuous state collocation nodes, collocation matrix and basis coefficients
s   = funnode(basis);      	% continuous state collocation nodes
Phi = funbase(basis,s);      % collocation matrix
c = zeros(nc,ni,nj);        % basis coefficients
for i=1:ni
for j=1:nj
  c(:,i,j) = Phi\v(:,i,j); 
end
end

% Solve finite horizon model by backward recursion.
if T<inf
  if output, display('Solving finite horizon model by backward recursion.'), end
  cc = zeros(nc,ni,T+2);
  vv = zeros(nc,ni,nj,T+2);
  xx = zeros(nc,dx,ni,nj,T+1);
  cc(:,:,T+2) = c;
  vv(:,:,:,T+2) = v;
  for t=T+1:-1:1
    [v,x] = vmax(model,basis,s,x,c);        % update optimal values and actions
    c = Phi\reshape(v,nc,ni*nj);            % update basis coefficients
    c = reshape(c,nc,ni,nj);                % reshape basis coefficient array
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

% Solve infinite-horizon model by function iteration.
if T==inf && strcmp(algorithm,'funcit')
  if output, display('Solving infinite-horizon model by function iteration.'), end
  for it=1:maxit
    cold = c;                                           % store old basis coefficients
    [v,x] = vmax(model,basis,s,x,c);                    % update optimal contingent values and actions
    c = Phi\reshape(v,nc,ni*nj);                        % update basis coefficients
    c = reshape(c,nc,ni,nj);                            % reshape basis coefficient array
    change = norm(c(:)-cold(:),inf);                    % compute change
    if output, fprintf ('%4i %10.1e\n',it,change), end	% print change
    if change<tol, break, end;                          % convergence check
  end
end

% Solve infinite horizon model by Newton's method.
if T==inf && strcmp(algorithm,'newton')
  display('Solving infinite horizon model by Newton''s method.')
  Phik = kron(eye(ni*nj),Phi);
  for it=1:maxit
    cold = c;                                           % store old basis coefficients
    [v,x,vc] = vmax(model,basis,s,x,c);                 % update optimal contingent values and actions
    c = c(:)-(Phik-vc)\(Phik*c(:)-v(:));                % update basis coefficients
    c = reshape(c,nc,ni,nj);                            % reshape basis coefficient array
    change = norm(c(:)-cold(:),inf);                    % compute change
    if output, fprintf ('%4i %10.1e\n',it,change), end	% print change
    if change<tol, break, end;                          % convergence check
  end
end

if output
  if it==maxit, display('Failure to converge in dpsolve2.'), end
  fprintf('Elapsed Time = %7.2f Seconds\n',toc)
end

% Check whether continuous state transitions remain within approximation interval
snmin =  inf;
snmax = -inf;
for i=1:ni
  [~,jmax] = max(v(:,i,:),[],3);
  for j=1:nj
    is = find(jmax==j);
    if isempty(is), continue, end
    ns = length(is);
    for k=1:size(e,1);
      ee = e(k+zeros(ns,1),:);
      sn = feval(func,'g',s(is,:),x(is,:,i,j),i,j,ee,params{:});
      snmin = min(snmin,min(sn));
      snmax = max(snmax,max(sn));
    end
  end
end
if output
  if any(snmin<basis.a-eps), display('Extrapolating below smin. Decrease smin.'), end;
  if any(snmax>basis.b+eps), display('Extrapolating above smax. Increase smax.'), end;
  disp([basis.a' snmin' snmax' basis.b'])
end

% Compute values and actions on refined continuous state array, if requested.
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

% Compute Bellman equation residual on refined continuous state array, if requested
if nargout>4
  vrproj = zeros(ns,ni,nj);
  for i=1:ni
    for j=1:nj
      vrproj(:,i,j) = funeval(c(:,i,j),basis,sr);
    end
  end
  resid = max(vrproj,[],3) - max(vr,[],3);
  resid = squeeze(resid);
end

% Eliminate singleton dimensions to facilitate analysis in calling program
c  = squeeze(c);
vr = squeeze(vr);
xr = squeeze(xr);

  
%% VMAX
%
%   Function called by dpsolve2 to maximize the optimand embedded in the
%   Bellman equation with respect to the continuous action x at ns input
%   continuous state nodes s, given a single discrete state i and a single
%   discrete action j.  Does so by solving associated Karush-Kuhn-Tucker
%   necessary conditions as a nonlinear complementarity problem, using a
%   derivative-based adapted Newton method (see Miranda and Fackler).
%   Function optionally also computes derivative of the optimand with
%   respect to the value function approximant basis function coefficients,
%   which is required to solve the collocation equation using Newton's
%   method.
%
% Usage
%   [v,x,vc] = vmax(model,basis,s,x,c)
% Let
%   n  = number of basis functions and continuous state collocation nodes
%   ns = number of continuous state nodes on refined output array
%   ds = dimension of continuous state s
%   dx = dimension of continuous action x
%   de = dimension of continuous state transition shock e
%   ni = number of discrete states i
%   nj = number of discrete actions j
%   nx = number of discretized continuous actions, if applicable
% Input
%   model  : structured array containing model specifications 
%   basis  : n-dimensional basis for real-valued functions defined on the 
%            continuous state space
%   s      : ns.ds array of continuous state nodes
%   x      : n.dx.ni.nj array of initial guesses for optimal continuous
%            actions at the ns continuous state nodes, per discrete state
%            and discrete action
%   c      : n.ni.nj array of discrete-action-contigent value function 
%            approximant basis function coefficients, per discrete state
%            and action
% Output
%   v      : n.ni.nj array of discrete-action-contingent values at the ns
%            continuous state nodes, per discrete state and discrete action
%   x      : ns.dx.ni.nj array of optimal continuous actions at the ns
%            continuous state nodes, per discrete state and discrete action
%   vc     : nn.nn (nn=nc*ni*nj) array of derivatives of values with
%            respect to basis coefficients


function [v,x,vc] = vmax(model,basis,s,x,c)

% Unpack user options
tol       = optget('dpsolve2','tol',sqrt(eps));
maxitncp  = optget('dpsolve2','maxitncp',50);  
ncpmethod = optget('dpsolve2','ncpmethod','minmax');

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
        [xl,xu] = feval(func,'b',s,[],i,j,[],params{:});
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
          vv(is,ix) = feval(func,'f',s(is,:),xx(is,:),i,j,[],params{:});
          for k=1:length(w)
            % Compute states next period and derivatives
            ee = e(k+zeros(length(is),1),:);
            snext = feval(func,'g',s(is,:),xx(is,:),i,j,ee,params{:});
            snext = real(snext);
            vn = zeros(length(is),nj);
            for in=1:ni
              if q(i,in,j)==0, continue, end
              prob = w(k)*q(i,in,j);
              for jn=1:nj
                vn(:,jn) = funeval(c(:,in,jn),basis,snext);
              end
              vv(is,ix) = vv(is,ix) + delta*prob*max(vn,[],2);
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
      [xl,xu] = feval(func,'b',s,[],i,j,[],params{:});
      for it=1:maxitncp
        % Initialize optimand and derivatives with reward terms
        [vv,vx,vxx] = feval(func,'f',s,x(:,:,i,j),i,j,[],params{:});
        for k=1:length(w)
          % Compute states next period and derivatives
          ee = e(k+zeros(ns,1),:);
          [snext,snx,snxx] = feval(func,'g',s,x(:,:,i,j),i,j,ee,params{:});
          snext = real(snext);
          B = funbasex(basis,snext,[0;1;2]);
          vn = zeros(ns,nj);
          vns = zeros(ns,nj,ds);
          vnss = zeros(ns,nj,ds,ds);
          for in=1:ni
            if q(i,in,j)==0, continue, end
            prob = w(k)*q(i,in,j);
            for jn=1:nj
              vn(:,jn)    = funeval(c(:,in,jn),basis,B);
              vns(:,jn,:) = funeval(c(:,in,jn),basis,B,eye(ds));
              for is=1:ds
                for js=is:ds
                  order = zeros(1,ds);
                  order(is) = order(is)+1;
                  order(js) = order(js)+1;
                  vnss(:,jn,is,js) = funeval(c(:,in,jn),basis,B,order);
                  vnss(:,jn,js,is) = vnss(:,jn,is,js);
                end
              end
            end
            [vn,jmax] = max(vn,[],2);
            vv = vv + delta*prob*vn;
            for jn=2:nj
              ins = find(jmax==jn);
              vns(ins,1,:)    = vns(ins,jn,:);
              vnss(ins,1,:,:) = vnss(ins,jn,:,:);
            end
            vns = reshape(vns(:,1,:),ns,ds);
            vnss = reshape(vnss(:,1,:,:),ns,ds,ds);
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
  vc = zeros(ns,ni,nj,ns,ni,nj);
  for i=1:ni
    for j=1:nj
      for k=1:length(w)
        ee = e(k+zeros(ns,1),:);
        snext = feval(func,'g',s,x(:,:,i,j),i,j,ee,params{:});
        B = funbase(basis,snext);
        vn  = zeros(ns,ni,nj);
        for in=1:ni
          for jn=1:nj
            vn(:,in,jn) = B*squeeze(c(:,in,jn));
          end
        end
        [vn,jmax] = max(vn,[],3);
        for in=1:ni
          if q(i,in,j)==0, continue, end
          prob = w(k)*q(i,in,j);
          vnc = zeros(ns,ns,nj);
          for jn=1:nj
            is = find(jmax(:,in)==jn);
            if isempty(is), continue, end
            vnc(is,:,jn) = B(is,:);
          end
          vnc = reshape(vnc,ns,1,1,ns,1,nj);
          vc(:,i,j,:,in,:) = vc(:,i,j,:,in,:) + delta*prob*vnc;
        end
      end
    end
  end
  vc = reshape(vc,ni*nj*ns,ni*nj*ns);
else
  vc = zeros(ns,ns);
  for k=1:length(w)
    ee = e(k+zeros(ns,1),:);
    snext = feval(func,'g',s,x,1,1,ee,params{:});
    vc = vc + delta*w(k)*funbase(basis,snext);
  end
end