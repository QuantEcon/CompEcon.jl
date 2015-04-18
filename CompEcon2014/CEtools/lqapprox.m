%% LQAPPROX
% 
%  Solves linear-quadratic approximation of discrete-time continuous-state 
%  continuous-action infinite-horizon single-agent dynamic optimization 
%  model solved by dpsolve. Uses QZ decomposition to solve Ricatti equation. 
%
%  Usage
%    [vlq,xlq,plq,sstar,xstar,pstar] = lqapprox(model,s,s0,x0)
%  Let
%    ns = number of continuous state nodes input
%    ds = dimension of continuous state s
%    dx = dimension of continuous action x
%  Input
%    model  : model structure
%    s      : ns.ds state nodes
%    s0     : 1.ds Taylor expansion state
%    x0     : 1.dx Taylor expansion action 
%  Output
%    vlq    : ns.1  optimal lq approximation values at state nodes
%    xlq    : ns.dx optimal lq approximation actions at state nodes
%    plq    : ns.ds optimal lq approximation shadow prices at state nodes
%    sstar  : 1.ds  computed steady-state state for lq approximation
%    xstar  : 1.dx  computed steady-state action for lq approximation
%    pstar  : 1.ds  computed steady-state shadow price for lq approximation
% 
%  Function File
%    A user-supplied function that returns the bound, reward, and state 
%    transition functions and their first and second derivatives with 
%    respect to the continuous action x, at an arbitrary number ns 
%    of continuous states s and continuous actions x, for a given discrete 
%    state i and a given discrete action j, according to the format
%      [out1,out2,out3] = func(flag,s,x,i,j,e,params)
%    Function File Input
%      flag      : flag that specifies the desired function
%      s         : ns.ds array of continuous state nodes
%      x         : ns.dx array of continuous actions
%      i         : ns.1 or 1.1 discrete state indices, integers between 1-ni
%      j         : a single discrete action index, an integer between 1-nj
%      e         : ns.de array of continuous state transition shocks
%      in        : ns.1 or 1.1 discrete state indices, integers between 1-ni
%      params    : user-supplied list of function parameters
%    Function File Output
%      if flag = 'b', returns lower and upper bounds on continuous action x
%        out1    : ns.dx lower bounds on continuous action x
%        out2    : ns.dx upper bounds on continuous action x
%        out3    : empty
%      if flag = 'f', returns reward and derivatives
%        out1    : ns.1 reward f
%        out2    : ns.dx first derivative of f with respect to x
%        out3    : ns.dx.dx second derivative of f with respect to x
%      if flag = 'g', returns continuous state transition and derivatives
%        out1    : ns.ds continuous state transtion g
%        out2    : ns.ds.dx first derivative of g with respect to x
%        out3    : ns.ds.dx.dx second derivative of g with respect to x
%      if flag = 'h', returns discrete state transition
%        out1    : ns.1 discrete state transtion h
%        out2    : empty
%        out3    : empty
%  See
%    dpsolve

%  Copyright(c) 1997-2010
%    Mario J. Miranda - miranda.4@osu.edu
%    Paul L. Fackler  - paul_fackler@ncsu.edu

function [vlq,xlq,plq,sstar,xstar,pstar] = lqapprox(model,s,s0,x0)

% Set model structure parameters to defaults, if nonexistent
if ~isfield(model,'ds'), model.ds = 1; end   % dimension of state s
if ~isfield(model,'dx'), model.dx = 1; end   % dimension of action x
if ~isfield(model,'e'),  model.e  = 0; end   % shocks
if ~isfield(model,'w'),  model.w  = 1; end   % shock probabilities
  
% Unpack model structure
func   = model.func;        % name of function file
params = model.params;      % parameters passed to function file
delta  = model.discount;    % discount factor
ds     = model.ds;          % dimension of continuous state s
dx     = model.dx;          % dimension of continuous action x
e      = model.e;           % shocks
w      = model.w;           % shock probabilties

ns = size(s,1);

% Fix shock at mean
estar = w'*e;

% Get derivatives
[f0,fx,fxx] = feval(func,'f',s0,x0,1,1,1,estar,params{:});
[g0,gx]     = feval(func,'g',s0,x0,1,1,1,estar,params{:});
fs  = fjac(func,2,'f',s0,x0,1,1,1,estar,params{:});
fxs = fjac(func,[2,2],'f',s0,x0,1,1,1,estar,params{:});
fss = fhess(func,2,'f',s0,x0,1,1,1,estar,params{:});
gs  = fjac(func,2,'g',s0,x0,1,1,1,estar,params{:});

% Reshape to ensure conformability
s0 = s0';
x0 = x0';
fs  = reshape(fs , 1,ds);
fx  = reshape(fx , 1,dx);
fss = reshape(fss,ds,ds);
fxs = reshape(fxs,dx,ds);
fxx = reshape(fxx,dx,dx);
g0  = reshape(g0 ,ds, 1);
gx  = reshape(gx ,ds,dx);
gs  = reshape(gs ,ds,ds);
fsx = fxs';
f0 = f0 - fs*s0 - fx*x0 + 0.5*s0'*fss*s0 + s0'*fsx*x0 + 0.5*x0'*fxx*x0;
fs = fs - s0'*fss - x0'*fxs;
fx = fx - s0'*fsx - x0'*fxx;
g0 = g0 - gs*s0 - gx*x0;

% Solve Ricatti equation using QZ decomposition
A = [eye(ds)          zeros(ds,dx+ds);
     zeros(dx,ds+dx) -delta*gx'      ;
     zeros(ds,ds+dx)  delta*gs'     ];
B = [ gs   gx  zeros(ds,ds);
     fsx' fxx  zeros(dx,ds);
    -fss -fsx  eye(ds)     ];
[S,T,Q,Z] = qzordered(A,B);
C  = real(Z(ds+1:end,1:ds)/Z(1:ds,1:ds));
X  = C(1:dx,:);
P  = C(dx+1:end,:);

% Compute steady-state state, action, and shadow price
t = [fsx' fxx delta*gx';fss fsx delta*gs'-eye(ds);gs-eye(ds) gx zeros(ds,ds)]\[-fx';-fs'; -g0];
sstar = t(1:ds);
xstar = t(ds+1:ds+dx);
pstar = t(ds+dx+1:ds+dx+ds);
vstar = (f0+fs*sstar+ fx*xstar+0.5*sstar'*fss*sstar+sstar'*fsx*xstar+0.5*xstar'*fxx*xstar)/(1-delta);

% Compute lq-approximation optimal policy and shadow price functions at
% state nodes
sstar = sstar';
xstar = xstar';
pstar = pstar';
s = s-sstar(ones(ns,1),:);
xlq = xstar(ones(1,ns),:) + s*X';
plq = pstar(ones(1,ns),:) + s*P';
vlq = vstar + s*pstar' + 0.5*sum(s.*(s*P'),2); 