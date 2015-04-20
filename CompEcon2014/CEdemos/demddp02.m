function demddp02

%% DEMDDP02 Asset replacement model

% Preliminary tasks
demosetup(mfilename)


%% Formulation

% Model Parameters
maxage  = 5;                  % maximum machine age
repcost = 75;               	% replacement cost
delta   = 0.9;                % discount factor

% State Space
S = (1:maxage)';             	% machine age
n = length(S);              	% number of states

% Action Space (keep=1, replace=2)
X = (1:2)';                  	% vector of actions
m = length(X);               	% number of actions

% Reward Function
f = [50-2.5*S-2.5*S.^2 (50-repcost)*ones(n,1)];
f(end,1) = -inf;

% State Transition Function
g = zeros(n,m);
for i=1:n
  g(i,1) = min(i+1,n);      	% keep
  g(i,2) = 1;               	% replace
end

% Model Structure
model.reward     = f;
model.transfunc  = g;
model.discount   = delta;


%% Solution

% Solve Bellman Equation
[v,x,pstar] = ddpsolve(model);
   

%% Analysis

% Plot Optimal Value
figure
plot(S,v)
title('Optimal Value Function')
xlabel('Age of Machine')
ylabel('Value')

% Simulate Model
sinit = min(S); nyrs = 12;
spath = ddpsimul(pstar,sinit,nyrs);

% Plot State Path
figure
plot(0:nyrs,S(spath))
title('Optimal State Path')
xlabel('Year')
ylabel('Age of Machine')
xlim([0 12])

% Save Plots as EPS Files
printfigures(mfilename,2)