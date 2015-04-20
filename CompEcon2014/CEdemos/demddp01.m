function demddp01

%% DEMDDP01 Mine management model

% Preliminary tasks
demosetup(mfilename)


%% Formulation

% Model Parameters
price = 1;                    % price of ore
sbar  = 100;                  % initial ore stock
delta = 0.9;                  % discount factor

% State Space
S = (0:sbar)';                % vector of states
n = length(S);                % number of states

% Action Space
X = (0:sbar)';                % vector of actions
m = length(X);                % number of actions

% Reward Function
f = zeros(n,m);
for i=1:n
  for k=1:m
    if X(k)<=S(i)
      f(i,k) = price*X(k)-(X(k)^2)./(1+S(i));
    else
      f(i,k) = -inf;
    end
  end
end

% State Transition Function
g = zeros(n,m);
for i=1:n
  for k=1:m
    snext = S(i)-X(k);
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
title('Optimal Extraction Policy')
xlabel('Stock')
ylabel('Extraction')

% Plot Value Function
figure
plot(S,v)
title('Optimal Value Function')
xlabel('Stock')
ylabel('Value')

% Simulate Model
sinit = max(S); nyrs = 15;
spath = ddpsimul(pstar,sinit,nyrs);

% Plot State Path
figure
plot(0:nyrs,S(spath))
title('Optimal State Path')
xlabel('Year')
ylabel('Stock')

% Save Plots as EPS Files
printfigures(mfilename,3)