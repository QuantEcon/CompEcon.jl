%% DPSTSTS
%
%  Computes invariant distribution of stochastic discrete-time 
%  continuous-state continuous-action single-agent dynamic optimization 
%  model solved by dpsolve.
%
%  Usage
%    [pp,ss,xx] = dpststs(model,nsmooth,nbin,s,x)
%  Let
%    ns = number of continuous state nodes on refined output grid
%    ds = dimension of continuous state s
%    dx = dimension of continuous action x
%  Input
%    model   : model structure
%    nsmooth : positive integer, smoothing parameter
%    nbin    : positive integer, number of nodes in distribution approximation
%    s       : ns.ds refined continuous state nodes
%    x       : ns.dx optimal continuous actions at refined continuous 
%                    state nodes, per discrete state and action
%  Output
%    sx      : midpoints of distribution histogram
%    pp      : associated invariant probabilities
%    xx      : control values at smid
% 
%  Note
%    Implemented only for 1-dimensional state space.
%  See
%    dpsolve

%  Copyright(c) 1997-2010
%    Mario J. Miranda - miranda.4@osu.edu
%    Paul L. Fackler  - paul_fackler@ncsu.edu

function [pp,ss,xx] = dpststs(model,nsmooth,nbin,s,x)

% Set shock distribution, number of discrete states and actions
% to defaults, if nonexistent or not passed as input
if ~isfield(model,'e'),  model.e = 0; end
if ~isfield(model,'w'),  model.w = 1; end
if ~isfield(model,'ds'), model.ds = 1; end  
if ~isfield(model,'dx'), model.dx = 1; end  
if ~isfield(model,'ni'), model.ni = 1; end
if ~isfield(model,'nj'), model.nj = 1; end
if ~isfield(model,'dx'), model.dx = 1; end

% Unpack model structure
func   = model.func;        % name of function file
params = model.params;      % parameters passed to function file
ds     = model.ds;          % dimension of continuous state s
dx     = model.dx;          % dimension of continuous action x
ni     = model.ni;          % number of discrete states
nj     = model.nj;          % number of discrete actions
e      = model.e;           % shocks
w      = model.w;           % shock probabilities

if ds>1||ni>1||nj>1||dx>1
  error('Implemented only for 1-dimensional continous state/action models'); 
end

smin = s(1);
smax = s(end);
swid = (smax-smin)./nbin;
smid = linspace(smin+swid/2,smax-swid/2,nbin)';

m  = length(w);
nobs = m^nsmooth;
ss = gridmake(smid);
pp = ones(size(ss));
for t=1:nsmooth
  xx = interp1(s,x,ss);
  [ss,xx,ee] = gridmake([ss xx],e);
  ss = feval(func,'g',ss,xx,1,1,ee,params{:});
  pp = pp*w'; 
  pp = pp(:);
end
ss = reshape(ss,nbin,nobs);
pp  = reshape(pp,nbin,nobs);

pi = zeros(nbin,nbin);
for j=1:nobs
  i = ceil((ss(:,j)-smin)/swid);
  i = min(max(i,1),nbin);
  ind=(1:nbin)'+(i-1)*nbin;
  pi(ind) = pi(ind) +  pp(i,j);
end

ss = gridmake(smid);
pp = markov(pi);
for t=1:nsmooth
  xx = interp1(s,x,ss);
  [ss,xx,ee] = gridmake([ss xx],e);
  ss = feval(func,'g',ss,xx,1,1,ee,params{:});
  pp = pp*w'; 
  pp = pp(:);
end