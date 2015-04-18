%% DDPSOLVE  
%
%  Discrete-Time Discrete-State Discrete-Action Single-Agent Dynamic 
%  Optimization Model Bellman Equation Solver
%
%  Computes solution to the Bellman equation of a discrete-time dynamic 
%  optimization model with discrete state i and discrete action x: 
%     V(i) = max_x {f(i,x)+delta*sum_j P_ij(x)V(j)}.
%  Returns optimal values v and optimal actions x.  
%
%  Usage
%    [v,x,pstar] = ddpsolve(model,v)
%  Let
%    n = number of discrete states
%    m = number of discrete actions
%    T = number of decision periods, if finite horizon
%  Input
%    model   : model structure variable
%    v       : initial guess for value function
%  Output
%    v       : optimal values (n.1 or n.T+1)
%    x       : optimal actions (n.1 or n.T)
%    pstar   : optimal transition probability matrix (n.n or n.n.T+1)
%
%  Model Structure Fields
%    discount  : discount factor
%    reward    : n.m matrix of reward values
%    transfunc : n.m transition function for deterministic problem 
%    transprob : m.n.n transition probabilities for stochastic problem
%    T         : number of decision periods, if finite horizon
%    vterm     : n.1 terminal values, if finite horizon
%  Note: Either transfunc or transprob, but not both, must be specified.
%  if T is omitted, an infinite horizon problem is assumed; if vterm is 
%  omitted, the terminal value if assumed to be 0.
%
%  Options
%    These options may be set by user with OPTSET (defaults in parentheses):
%    algorithm : solution algorithm, either Newton method ('newton') or 
%                function iteration 'funcit'
%    tol       : solution algorithm convergence tolerance (sqrt(eps))
%    maxit     : solution algorithm maximum number of iterations (500)

%  Copyright(c) 1997-2010
%    Mario J. Miranda - miranda.4@osu.edu
%    Paul L. Fackler  - paul_fackler@ncsu.edu

function [v,x,pstar] = ddpsolve(model,v)

% Set options to defaults, if not set by user with OPTSET (see above)
maxit = optget('ddpsolve','maxit',100);
tol   = optget('ddpsolve','tol',sqrt(eps));

% Unpack model structure
delta = model.discount;
f     = model.reward;
if isfield(model,'horizon')
  T = model.horizon;
else
  T = inf;
end
if isfield(model,'transfunc')
  P = expandg(model.transfunc);
else
  P = model.transprob;
  [n,m] = size(f);
  if length(size(P)) == 3 && all(size(P) == [m n n])
    P = permute(P,[2 1 3]);
    P = reshape(P,m*n,n);
  else
    error('P has the wrong size')
  end
end

if T <inf
  algorithm = 'back';
else
  algorithm = optget('ddpsolve','algorithm','newton');
end

n = size(f,1);
if T<inf
  if isfield(model,'vterm')
    v = model.vterm;
    if size(v,1) ~= n
      error('model.vterm has improper size')
    end
  else v = zeros(n,1);
  end
else
  if nargin<3, v = zeros(n,1); end
end

% Solve Bellman equation
switch algorithm
  case 'newton'
    disp('Solve Bellman equation via Newton method')
    for it=1:maxit
      vold = v;                             % store old value
      [v,x] = valmax(v,f,P,delta);          % update policy
      [pstar,fstar] = valpol(x,f,P);        % induced P and f
      v = (speye(n)-delta*pstar)\fstar;     % update value
      change = norm(v-vold);                % compute change
      fprintf ('%5i %10.1e\n',it,change)    % print progress
      if change<tol, break, end;            % convergence check
    end
    if change>tol, warning('Failure to converge in ddpsolve'), end;
  case 'funcit'
    disp('Solve Bellman equation via function iteration')
    for it=1:maxit
      vold = v;                             % store old value
      [v,x] = valmax(v,f,P,delta);          % update policy
      change = norm(v-vold);                % compute change
      fprintf ('%5i %10.1e\n',it,change)    % print progress
      if change<tol, break, end;            % convergence check
    end
    pstar = valpol(x,f,P);
    if change>tol, warning('Failure to converge in ddpsolve'), end;
  case 'back'
    disp('Solve Bellman equation via backward recursion')
    x = zeros(n,T);
    v = [zeros(n,T) v];
    pstar = zeros(n,n,T);
    for t=T:-1:1                                   % backward recursion
      [v(:,t),x(:,t)] = valmax(v(:,t+1),f,P,delta);  % Bellman equation
      pstar(:,:,t) = valpol(x(:,t),f,P);
    end
  otherwise
    error('algorithm must be newton or funcit')
end


function [v,x] = valmax(v,f,P,delta)
[n,m] = size(f);
[v,x] = max(f+delta*reshape(P*v,n,m),[],2);

function [pstar,fstar] = valpol(x,f,P)
n = size(f,1);
i = n*(x-1)+(1:n)';
fstar = f(i);
pstar = P(i,:);

function P = expandg(g)
[n,m] = size(g);
P = sparse(1:n*m,g(:),1,n*m,n);