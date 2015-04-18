%% DDPSIMUL  
%
%  Discrete-Time Discrete-State Discrete-Action Single-Agent Dynamic 
%  Optimization Model Simulator
%
%  Simulates the Markov chain followed by the optimized state variable of
%  a discrete-time discrete-state discrete-action single-agent dynamic 
%  optimization model solved using ddpsolve.
%
%  Usage
%    [spath,xpath] = ddpsimul(pstar,s,nper,x)
%  Let
%    n    = number of discrete states
%    nrep = number of replications
%    nper = number of periods simulated
%  Input
%    pstar  : n.1    optimal state transition matrix
%    s      : nrep.1 initial states
%    nper   : number of periods simulated, a positive integer
%    x      : n.1    optimal actions
%  Output
%    spath  : nrep.nper+1 simulated states
%    xpath  : nrep.nper+1 simulated actions
%  See
%    ddpsolve

%  Copyright(c) 1997-2010
%    Mario J. Miranda - miranda.4@osu.edu
%    Paul L. Fackler  - paul_fackler@ncsu.edu

function [spath,xpath] = ddpsimul(pstar,s,nper,x)

l = length(size(pstar));
n = size(pstar,2);
u = ones(n,1);
nrep = length(s);

if l == 2 % Infinite Horizon Model
  spath = zeros(nrep,nper+1);
  cp = cumsum(pstar,2);
  for t = 1:nper+1
    spath(:,t) = s;
    if t<= nper
      r = rand(nrep,1);
      s = 1+sum(r(:,u)>cp(s,:),2);
    end
  end
else    % Finite Horizon Model
  T = size(pstar,3);
  if nper>T,
    warning('Request for simulations beyond time horizon ignored')
  end
  nper = min(nper,T);
  spath = zeros(nrep,nper+1);
  for t = 1:nper+1
    spath(:,t) = s;
    if t<= nper
      cp = cumsum(pstar(:,:,t),2);
      r = rand(nrep,1);
      s = 1+sum(r(:,u)>cp(s,:),2);
    end
  end
end

if nargout>1
  xpath = zeros(nrep,nper+1);
  xpath(:) = x(spath(:));
end