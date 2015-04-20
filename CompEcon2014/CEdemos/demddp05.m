function demddp05

%% DEMDDP05 Water management model

% Preliminary tasks
demosetup(mfilename)


%% Formulation

% Model Parameters
alpha1 = 14;                  % producer benefit function parameter
beta1  = 0.8;                 % producer benefit function parameter
alpha2 = 10;                  % recreational user benefit function parameter
beta2  = 0.4;                 % recreational user benefit function parameter
maxcap = 30;                  % maximum dam capacity
r = [0 1 2 3 4];              % rain levels
p = [0.1 0.2 0.4 0.2 0.1];    % rain probabilities
delta  = 0.9;                 % discount factor

% State Space
S = (0:maxcap)';              % vector of states
n = length(S);                % number of states

% Action Space
X = (0:maxcap)';              % vector of actions
m = length(X);                % number of actions

% Reward Function
f = zeros(n,m);
for i=1:n
  for k=1:m
    if k>i
      f(i,k) = -inf;
    else
      f(i,k) = alpha1*X(k).^beta1+alpha2*(S(i)-X(k)).^beta2;
    end
  end
end

% State Transition Probability Matrix
P = zeros(m,n,n);
for k=1:m
  for i=1:n
    for j=1:length(r)
      snext = min(S(i)-X(k)+r(j),maxcap);
      inext = getindex(snext,S);
      P(k,i,inext) = P(k,i,inext) + p(j);
    end
  end
end

% Model Structure
model.reward     = f;
model.transprob  = P;
model.discount   = delta;


%% Solution

% Solve Bellman Equation
[v,x,pstar] = ddpsolve(model);
   

%% Analysis

% Plot Optimal Policy
figure
h=plot(S,X(x),'*'); % set(h,'FaceColor',[.75 .75 .75])
axis([0 maxcap -inf inf])
title('Optimal Irrigation Policy')
xlabel('Water Level')
ylabel('Irrigation')
xlim([-1 31])
ylim([0 6])

% Plot Value Function
figure
plot(S,v)
title('Optimal Value Function')
xlabel('Water Level')
ylabel('Value')

% Simulate Model
sinit = ones(10000,1);
nyrs  = 30;
spath = ddpsimul(pstar,sinit,nyrs);

% Plot State Path
figure
plot(0:nyrs,mean(S(spath)))
title('Optimal State Path')
xlabel('Year')
ylabel('Water Level')

% Compute Steady-State Distribution of Water Level
pi = markov(pstar);
figure
h=bar(S,pi,1); set(h,'FaceColor',[.75 .75 .75])
title('Steady State Distribution')
xlabel('Water Level')
ylabel('Probability')
xlim([-1 31])
ylim([0 0.16])

% Compute Steady-State Mean Water Level
avgstock = pi'*S;
fprintf('\nSteady-state Stock        %8.2f\n',avgstock)

% Save Plots as EPS Files
printfigures(mfilename,4)