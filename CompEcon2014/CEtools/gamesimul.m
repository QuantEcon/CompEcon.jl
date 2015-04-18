%% GAMESIMUL
%
%  Discrete-Time Continuous-State Continuous-Action Infinite-Horizon
%  Multi-Agent Dynamic Game Model Simulator
%
%  Simulates the Markov process followed by the optimized state variable of
%  a discrete-time continuous-state continuous-action infinite-horizon
%  multi-agent dynamic game model solved using gamesolve.
%
%  Usage
%    [ssim,xsim] = gamesimul(model,ss,nper,sres,xres)
%  Input
%    model   : model structure variable
%    ss      : k by d vector of initial states
%    nper    : number of simulated time periods
%    sres    : coordinates of the evaluation grid (from dpsolve)
%    xres    : optimal control function values at grid defined by sres
%  Output
%    ssim    : k by d by nper+1 vector of simulated states       
%    xsim    : k by d by nper+1 vector of simulated actions
%    For finite horizon problems, xsim is k by d by nper
%    If d=1, ssim and xsim are k by nper+1
%  Uses
%    discrand
%  See
%    gamesolve

%  Copyright(c) 1997-2010
%    Mario J. Miranda - miranda.4@osu.edu
%    Paul L. Fackler  - paul_fackler@ncsu.edu

function [ssim,xsim] = gamesimul(model,ss,nper,sres,xres)

% Dummy shocks and weidhts if model is deterministic
if isfield(model,'e'), e=model.e; else e=0; end;
if isfield(model,'w'), w=model.w; else w=1; end;

func   = model.func;
params = model.params;

nrep = size(ss,1);
ds   = size(ss,2);
st   = gridmake(sres);
dx   = ds*length(xres(:))/length(st(:));
ssim = zeros(nrep,ds,nper+1);
xsim = zeros(nrep,dx,nper+1);

nx = numel(xres)/dx;
xres = reshape(xres,nx,dx);
for t=1:nper+1
  xx = minterp(sres,xres,ss);
  ssim(:,:,t) = ss;
  xsim(:,:,t) = xx;
  ee = e(discrand(nrep,w),:);
  ss = feval(func,'g',1,ss,xx,ee,params{:});
end

if ds==1; ssim=squeeze(ssim); end
if dx==1; xsim=squeeze(xsim); end