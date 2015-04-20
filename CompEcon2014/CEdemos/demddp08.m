function demddp08

%% DEMDDP08 Job search model

% Preliminary tasks
demosetup(mfilename)


%% Formulation

% Model Parameters
u     =  50;                  % weekly unemp. benefit
v     =  60;                  % weekly value of leisure
pfind = 0.90;                 % prob. of finding job
pfire = 0.10;                 % prob. of being fired
delta = 0.99;                 % discount factor

% State Space
S = (1:2)';                   % vector of states
n = length(S);                % number of states

% Action Space (idle=1, active=2)
X = (1:2)';                  	% vector of actions
m = length(X);               	% number of actions

% Reward Function
f = zeros(n,m);
f(:,1) = v;                   % gets leisure
f(1,2) = u;                   % gets benefit

% State Transition Probability Matrix
P = zeros(n,n,n);
P(1,:,1) = 1;                 % remains unemployed
P(2,1,1) = 1-pfind;           % finds no job
P(2,1,2) = pfind;             % finds job
P(2,2,1) = pfire;             % gets fired
P(2,2,2) = 1-pfire;           % keeps job

% Model Structure
model.reward     = f;
model.transprob  = P;
model.discount   = delta;


%% Solution

% Solve Bellman Equation
xtable = [];
wage   = 55:65;
for w=wage
  f(2,2) = w; model.reward = f;  % vary wage
  [v,x]  = ddpsolve(model);      % solve via policy iteration
  xtable = [xtable x];           % tabulate
end


%% Analysis

% Display Optimal Policy
fprintf('\nOptimal Job Search Strategy')
fprintf('\n  (1=innactive, 2=active)\n')
fprintf('\nWage  Unemployed  Employed\n')
fprintf('%4i  %10i%10i\n',[wage;xtable])