function demddp09

%% DEMDDP09 Deterministic cow replacement model

% Preliminary tasks
demosetup(mfilename)


%% Formulation

% Model Parameters
delta = 0.9;                  % discount factor
cost  = 500;                  % replacement cost
price = 150;                  % milk price

% State Space
S = (1:10)';                  % lactation states
n = length(S);                % number of states

% Action Space (keep='K', replace='R')
X = ['K';'R'];                % keep or replace
m = length(X);                % number of actions

% Reward Function
f = zeros(n,m);
y = (-0.2*S.^2+2*S+8);        % yield per lactation
f = [price*y price*y-cost];   % net revenue by action
f(10,1) = -inf;               % force replace at lactation 10

% State Transition Function
g = zeros(n,m);
for i=1:n
  g(i,1) = min(i+1,n);       % Raise lactation number by 1, if keep
  g(i,2) = 1;                % Lactation number reverts to 1, if replace
end

% Model Structure
model.reward     = f;
model.transfunc  = g;
model.discount   = delta;


%% Solution

% Solve Bellman Equation
[v,x] = ddpsolve(model);
   

%% Analysis

% Plot Optimal Policy
figure
bar(S,x)
xlabel('Age')
ylabel('Optimal Decision')
set(gca,'YTick',[1 2])
set(gca,'YTickLabel',{'Keep','Replace'})

% Plot Value Function
figure
plot(S,v)
xlabel('Age')
ylabel('Optimal Value')

% Save Plots as EPS Files
printfigures(mfilename,2)