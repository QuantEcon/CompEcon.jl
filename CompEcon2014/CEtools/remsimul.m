%% REMSIMUL
%
%  Discrete-Time Continuous-State Continuous-Response Infinite-Horizon 
%  Rational Expectations Model Simulator
%
%  Simulates the Markov process followed by the state variable of a 
%  discrete-time continuous-state continuous-response infinite-horizon
%  rational expectations model solved using remsolve.
%
% Usage
%   [ssim,xsim] = remsimul(model,basis,nper,sinit,s,v,x)
% Let
%   ns   = number of continuous state nodes
%   ds   = dimension of state variable s
%   dx   = dimension of response variable x
%   de   = dimension of shock variable e
%   nrep = number of replications
%   nper = number of periods simulated
% Input
%   model  : model structure
%   basis  : approximation structure
%   nper   : number of periods simulated, a positive integer
%   sinit  : nrep.ds array of nrep initial continuous states
%   v      : ns.1  optimal values at state collocation nodes
%   x      : ns.dx responses at state collocation nodes
% Output
%   ssim   : nrep.nper+1.ds array of simulated continuous states
%   xsim   : nrep.nper+1.dx array of simulated continuous actions
%   if nrep=1, arrays are squeezed along this singleton dimension
%  See
%    remsolve
%
% Copyright(c) 1997-2010
%   Mario J. Miranda - miranda.4@osu.edu
%   Paul L. Fackler  - paul_fackler@ncsu.edu

function [ssim,xsim] = remsimul(model,basis,nper,sinit,s,x)

% Set shock distribution at default, if nonexistent
if ~isfield(model,'e'),  model.e = 0; end
if ~isfield(model,'w'),  model.w = 1; end
if ~isfield(model,'ds'), model.ds = 1; end
if ~isfield(model,'dx'), model.dx = 1; end
  
% Unpack model structure
func   = model.func;        % name of function file
params = model.params;      % parameters passed to function file
ds     = model.ds;          % dimension of continuous state s
dx     = model.dx;          % dimension of continuous action x
e      = model.e;           % shocks
w      = model.w;           % shock probabilities

% Determine number of collocation nodes and variable dimensions
ns = size(s,1);             % number of continuous state variables input
nrep = size(sinit,1);       % number of replications

if nargin<6
  x = zeros(ns,dx);
end

% Reshape x and v to conform to format recognized by this routine
x = reshape(x,ns,dx);

% Compute collocation matrix and basis coefficients for optimal policy and 
% value approximations, per discrete state i and discrete action j
Phi = funbase(basis,s);
cx  = Phi\x;      

% Initialize output arrays
ssim = zeros(nrep,nper+1,ds);
xsim = zeros(nrep,nper+1,dx);
ss = sinit;
for ip=1:nper+1
  xx = funeval(cx,basis,ss);
  [xl,xu] = feval(func,'b',ss,xx,[],[],[],params{:});
  xx = min(max(xx,xl),xu);
  ssim(:,ip,:) = ss;
  xsim(:,ip,:) = xx;
  if ip<nper+1
    ee = e(discrand(nrep,w),:);
    ss = feval(func,'g',ss,xx,[],[],ee,model.params{:});
  end
end
% if nrep==1
%   ssim = reshape(ssim,nper+1,ds);
%   xsim = reshape(xsim,nper+1,dx);
% end