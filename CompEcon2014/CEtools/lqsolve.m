%% LQSOLVE 
%
%  General Infinite-Horizon Linear-Quadratic Dynamic Optimization Solver.
%
%  Uses QZ decomposition to solve the Ricatti equation of a deterministic
%  stationary infinite-horizon linear-quadratic dynamic optimization model
%
%    max_x f0 + fs*s +  fx*x + 0.5*s'*fss*s + s'*fsx*x +0.5*x'*fxx*x
%    s.t. s' = g0 + gs*s + gx*x
%
%    The optimal policy function is
%       x(s) = xstar + X*(s-sstar)
%    The shadow price function is
%       p(s) = pstar + P*(s-sstar)
%    The value function is
%       V(s) = vstar + pstar*(s-sstar) + 0.5*(s-sstar)'*P*(s-sstar)
%    The controlled state process is
%      snext = sstar + G*(s-sstar)
%
%  Usage
%    [X,P,G,sstar,xstar,pstar,vstar] = lqsolve(f0,fs,fx,fss,fsx,fxx,g0,gs,gx,delta)
%  Let
%    ds = dimension of state s, a row vector
%    dx = dimension of action x, a row vector
%  Input
%    f0     : 1.1   objective function parameter
%    fs     : 1.ds  objective function parameter
%    fx     : 1.dx  objective function parameter
%    fss    : ds.ds objective function parameter
%    fsx    : ds.dx objective function parameter
%    fxx    : dx.dx objective function parameter
%    g0     : ds.1  state transition function parameter
%    gs     : ds.ds state transition function parameter
%    gx     : ds.dx state transition function parameter
%    delta  : discount factor
%  Output
%    X      : dx.ds response of optimal action to state
%    P      : ds.ds response of shadow price to state
%    G      : ds.ds optimal state transition response
%    sstar  : ds.1 steady-state state
%    xstar  : ds.1 steady-state action
%    pstar  : ds.1 steady-state shadow price
%    vstar  : steady-state value
%
% Note: A stochastic LQ problem with additive, zero-mean state tranistion 
% shock e will have the same optimal policy and shadow price function as
% the deterministic problem with the shock ignored.  The value function of
% the stochastic problem will equal the value function of the deterministic 
% model plus the constant 0.5*delta*Ee'Pe/(1-delta), were E is the
% expectation operator.

%  Copyright(c) 1997-2013
%    Mario J. Miranda - miranda.4@osu.edu
%    Paul L. Fackler  - paul_fackler@ncsu.edu

function [X,P,G,sstar,xstar,pstar,vstar] = lqsolve(f0,fs,fx,fss,fsx,fxx,g0,gs,gx,delta)

% Determine dimensions of state and action variables
ds = length(fs);
dx = length(fx);

% Check conformability
if size(fs)  ~= [ 1 ds], error('error in lqsolve, dimension of fs' ), end;
if size(fx)  ~= [ 1 dx]; error('error in lqsolve, dimension of fx' ), end;
if size(fss) ~= [ds ds]; error('error in lqsolve, dimension of fss'), end;
if size(fsx) ~= [ds dx]; error('error in lqsolve, dimension of fss'), end;
if size(fxx) ~= [dx dx]; error('error in lqsolve, dimension of fxx'), end;
if size(g0)  ~= [ds  1]; error('error in lqsolve, dimension of g0' ), end;
if size(gs)  ~= [ds ds]; error('error in lqsolve, dimension of gs' ), end;
if size(gx)  ~= [ds dx]; error('error in lqsolve, dimension of gx' ), end;

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
G  = gs+gx*X;

% Compute steady-state state, action, and shadow price
t = [fsx' fxx delta*gx';fss fsx delta*gs'-eye(ds);gs-eye(ds) gx zeros(ds,ds)]\[-fx';-fs'; -g0];
sstar = t(1:ds);
xstar = t(ds+1:ds+dx);
pstar = t(ds+dx+1:ds+dx+ds);
vstar = (f0+fs*sstar+ fx*xstar+0.5*sstar'*fss*sstar+sstar'*fsx*xstar+0.5*xstar'*fxx*xstar)/(1-delta);

% Alternate computation of X, P, and G using Matlab dare function.  This
% method fails with unacceptable frequency because spectrum is too close 
% to the unit circle.
% [P,L,X] = dare(gs,gx,fss/delta,fxx/delta,fsx/delta,eye(ds)/sqrt(delta));
% X = -X;
% G = gs+gx*X;