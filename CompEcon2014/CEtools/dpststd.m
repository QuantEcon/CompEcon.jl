%% DPSTSTD
%
%  Computes steady-state state, action, and shadow price of deterministic 
%  discrete-time continuous-state continuous-action single-agent dynamic 
%  optimization model solved by dpsolve.
%
%  Usage
%    [sstar,xstar,pstar] = dpststd(model,s0,x0,p0)
%  Input
%    model  : model structure
%    s0     : 1.ds initial guess for steady-state state
%    x0     : 1.dx initial guess for steady-state action
%    p0     : 1.ds initial guess for steady-state shadow price
% Output
%    sstar  : 1.ds steady-state state
%    xstar  : 1.dx steady-state action
%    pstar  : 1.ds steady-state shadow price
%  See
%    dpsolve

%  Copyright(c) 1997-2010
%    Mario J. Miranda - miranda.4@osu.edu
%    Paul L. Fackler  - paul_fackler@ncsu.edu

function [sstar,xstar,pstar] = dpststd(model,s0,x0,p0)

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

% Fix shock at mean
e = w'*e;

% Use Broyden's method to solve steady-state equilibrium conditions
theta = [s0(:);x0(:);p0(:)];
theta = broyden(@dpssres,theta,func,params,e,delta,ds,dx);
sstar = theta(1:ds)';
xstar = theta(ds+1:ds+dx)';
pstar = theta(ds+dx+1:end)';


%% DPSSRES Residual function for computing steady-state
function r = dpssres(theta,func,params,e,delta,ds,dx)

% Unpack theta
s = theta(1:ds)';
x = theta(ds+1:ds+dx)';
p = theta(ds+dx+1:end)';

% Get derivatives
[f0,fx] = feval(func,'f',s,x,1,1,e,params{:});
[g0,gx] = feval(func,'g',s,x,1,1,e,params{:});
fs      = fjac(func,2,'f',s,x,1,1,e,params{:});
gs      = fjac(func,2,'g',s,x,1,1,e,params{:});

% Reshape to ensure conformability
fs = reshape(fs, 1,ds);
fx = reshape(fx, 1,dx);
gx = reshape(gx,ds,dx);
gs = reshape(gs,ds,ds);

% Compute residual
[xl,xu] = feval(func,'b',s,x,1,1,e,params{:});
rx = min(max(fx+delta*p*gx,xl-x),xu-x);
r = [s-g0 rx p-(fs+delta*p*gs)]';