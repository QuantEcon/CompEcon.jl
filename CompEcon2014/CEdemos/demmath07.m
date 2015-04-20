function demmath07

%% DEMDP07 Operations with Markov Chains

% Preliminary tasks
demosetup(mfilename)


%% FORMULATION
    
% Model Parameters
delta = 0.90;        % discount factor
gamma = 0.07;        % unemployment rate
eta   = 2.0;         % expected duration of unemployment
y     = [0.5;1.0];   % income per employment state

% Employment State Transition Probabilities 
q      = zeros(2,2);
q(1,1) = 1-1/eta;
q(2,1) = gamma*(1-q(1,1))/(1-gamma);
q(1,2) = 1-q(1,1);
q(2,2) = 1-q(2,1);
% q = [0.1 0.9; 0.1 0.9]

% Employment State Stationary Disbritution and Expected Visit Durations
p   = markov(q);
eta = 1./(1-diag(q));

% Print Employment State Transition Probabilities, Stationary Disbritution and Expected Visit Durations
fprintf('\n\n')
fprintf('Employment State Transition and Stationary Probabilities and Expected State Durations\n')
fprintf('           Unemployed   Employed   Stationary\n')
fprintf('Unemployed    %7.3f  %9.3f  %11.3f\n',q(1,:)',p(1))
fprintf('Employed      %7.3f  %9.3f  %11.3f\n',q(2,:)',p(2))
fprintf('Duration      %7.3f  %9.3f\n',eta(:))

% Compute present value of expected income per state
e = (eye(2,2)-delta*q)\y

% Compute ergodic mean income
Ey = p'*y;

% Compute ergodic standard deviation of income
Sy = sqrt(p'*(y.^2)-Ey^2);

% Compute ergodic income autocorrelation
Ay = (((p.*y)'*q*y)-Ey^2)/(p'*(y.^2)-Ey^2);

% Simulate
n = 1000000;
i = discrand(n,p);
j = markovsim(i,q);

% Compute simulated ergodic mean incomes
m = mean([y(i);y(j)]);

% Compute simulated ergodic standard deviation of income
s = std([y(i);y(j)]);

% Compute simulated ergodic income autocorrelation
a = (sum(y(i).*y(j))/n-m^2)/(s^2);

fprintf('\n\n')
fprintf('Employment State Transition and Stationary Probabilities and Expected State Durations\n')
fprintf('           True Simulated\n')
fprintf('Mean                %7.3f  %7.3f\n',Ey,m)
fprintf('Standard Deviation  %7.3f  %7.3f\n',Sy,s)
fprintf('Autocorrelation     %7.3f  %7.3f\n',Ay,a)
