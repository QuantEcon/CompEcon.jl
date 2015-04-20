function demddp10

%% DEMDDP10 Stochastic cow replacement model

% Preliminary tasks
demosetup(mfilename)


%% Formulation

% Model Parameters
delta = 0.9;                  % discount factor
cost  = 500;                  % replacement cost
price = 150;                  % milk price

% State Space
s1 = (1:10)';                 % lactation states
s2 = [0.8;1.0;1.2];           % productivity states
n1 = length(s1);              % number of lactation states
n2 = length(s2);              % number of productivity states
[S1,S2] = gridmake(s1,s2);    % combined state grid
n = n1*n2;                    % total number of states

% Action Space (keep='K', replace='R')
X = ['K','R'];                % keep or replace
m = length(X);               	% number of actions

% Reward Function
y = (-0.2*S1.^2+2*S1+8).*S2;  % yield per lactation
f = [price*y price*y-cost];   % net revenue by action
f(S1==10,1) = -inf;           % force replace at lactation 10

% State Transition Probability Matrix
P = zeros(2,n1,n2,n1,n2);
for i=1:n1
  for j=1:n2
    if i<10
      P(1,i,j,i+1,j) = 1;     % Raise lactation number by 1, if keep
    else
      P(1,i,j,1,1) = 0.2;     % Forced replacement after lactation 10
      P(1,i,j,1,2) = 0.6;
      P(1,i,j,1,3) = 0.2;
    end
    P(2,i,j,1,1) = 0.2;       % Optional replacement
    P(2,i,j,1,2) = 0.6;
    P(2,i,j,1,3) = 0.2;
  end
end
P = reshape(P,2,n,n);

% Model Structure
model.reward     = f;
model.transprob  = P;
model.discount   = delta;


%% Solution

% Solve Bellman Equation
[v,x,pstar] = ddpsolve(model);
   

%% Analysis

% Display Optimal Policy
disp('Optimal Policy')
disp('     Age       Lo       Med        Hi')
fprintf('%8i %8c  %8c  %8c\n',[s1 reshape(X(x),n1,n2)]')

% Plot Value Function
figure
plot(s1,reshape(v,n1,n2))
xlabel('Age')
ylabel('Optimal Value')
legend('Low','Med','Hi')

% Compute Steady-State distribution
pi = markov(pstar);

% Display Steady-State distribution
disp('          Invariant Distribution     ')
disp('     Age       Lo       Med        Hi')
fprintf('%8i %8.3f  %8.3f  %8.3f\n',[s1 reshape(pi,n1,n2)]')

% Compute Steady-State Mean Cow Age and Productivity
avgage = pi'*S1;
avgpri = pi'*S2;
fprintf('\nSteady-state Age          %8.2f\n',avgage)
fprintf('\nSteady-state Productivity %8.2f\n',avgpri)

% Save Plots as EPS Files
printfigures(mfilename,1)