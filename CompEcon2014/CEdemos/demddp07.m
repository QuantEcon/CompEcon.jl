function demddp07

%% DEMDDP07 Renewable resource model

% Preliminary tasks
demosetup(mfilename)


%% Formulation

% Model Parameters
delta =  0.9;                  % discount factor
alpha =  4.0;                  % growth function parameter
beta  =  1.0;                  % growth function parameter
gamma =  0.5;                  % demand function parameter
cost  =  0.2;                  % unit cost of harvest

% State Space
smin = 0;                      % minimum state
smax = 8;                      % maximum state
n    = 200;                    % number of states
S    = nodeunif(n,smin,smax);  % vector of states

% Action Space
xmin = 0;                      % minimum action
xmax = 6;                      % maximum action
m    = 100;                    % number of actions
X    = nodeunif(m,xmin,xmax);  % vector of actions

% Reward Function
f = zeros(n,m);
for k=1:m
  f(:,k) = (X(k).^(1-gamma))/(1-gamma)-cost*X(k);
  f(S<X(k),k) = -inf;
end

% State Transition Function
g = zeros(n,m);
for i=1:n
  for k=1:m
    snext  = alpha*(S(i)-X(k)) - 0.5*beta*(S(i)-X(k)).^2;
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
plot(S,X(x))
title('Optimal Harvest Policy')
xlabel('Stock')
ylabel('Harvest')

% Plot Value Function
figure
plot(S,v)
title('Optimal Value Function')
xlabel('Stock')
ylabel('Value')

% Simulate Model
nyrs = 20;
spath = ddpsimul(pstar,n,nyrs);

% Plot State Path
figure
plot(0:nyrs,S(spath))
title('Optimal State Path')
xlabel('Year')
ylabel('Stock')

% Plot Optimal Transition Function
figure
[ii,jj]=find(pstar);
plot(S(ii),S(jj),S,S,'--')
title('Optimal State Transition Function')
xlabel('S(t)')
ylabel('S(t+1)')

% Save Plots as EPS Files
printfigures(mfilename,4)