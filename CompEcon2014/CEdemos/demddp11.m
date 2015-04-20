function demddp11

%% DEMDDP11 Stochastic optimal growth model

% Preliminary tasks
demosetup(mfilename)


%% Formulation

% Model Parameters
delta =  0.9;                 % discount factor
alpha =  0.2;                 % utility parameter
beta  =  0.5;                 % production parameter
gamma =  0.9;                 % capital survival rate

% State Space
smin  =  1.0;                 % minimum state
smax  = 10.0;                 % maximum state
n = 200;                      % number of states
S = nodeunif(n,smin,smax);    % vector of states

% Action Space
xmin  = 0.65;                 % minimum action
xmax  = 0.99;                 % maximum action
m =  500;                     % number of actions
X = nodeunif(m,xmin,xmax);    % vector of actions

% Reward Function
f = zeros(n,m);
for k=1:m
  f(:,k) = ((S-X(k)*S).^(1-alpha))/(1-alpha);
end

% State Transition Function
g = zeros(n,m);
for i=1:n
  for k=1:m
    snext = gamma*X(k)*S(i) + (X(k)*S(i)).^beta;
    g(i,k) = getindex(snext,S);
  end
end

% Model Structure
model.reward     = f;
model.transfunc  = g;
model.discount   = delta;


%% Solution

% Solve Bellman Equation
[v,x,pstar] = ddpsolve(model);
   

%% Analysis

% Plot Optimal Policy
figure
plot(S,X(x).*S)
title('Optimal Investment')
xlabel('Wealth')
ylabel('Investment')

% Plot Optimal Policy
figure
plot(S,S-X(x).*S)
title('Optimal Consumption')
xlabel('Wealth')
ylabel('Consumption')

% Plot Value Function
figure
plot(S,v)
title('Optimal Value Function')
xlabel('Wealth')
ylabel('Value')

% Simulate Model
nyrs = 20;
st = ddpsimul(pstar,smin,nyrs);

% Plot State Path
figure
plot(0:nyrs,S(st))
title('Optimal State Path')
xlabel('Year')
ylabel('Wealth')

% Compute Steady State Distribution and Mean
pi = markov(pstar);
avgstock = pi'*S;
fprintf('\nSteady-state Wealth     %8.2f\n',avgstock)

% Save Plots as EPS Files
printfigures(mfilename,4)