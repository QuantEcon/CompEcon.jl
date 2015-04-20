function demmath06

%% Simulate Simple Markov Chain

% Preliminary tasks
demosetup(mfilename)


%% FORMULATION
    
% Model Parameters
gamma    = 0.07;                % aggregate unemployment rate
eta      = 2.0;                 % expected duration of unemployment
y        = [0.5;1.0];           % income per employment state
delta    = 0.90;                % discount factor

% Employment Transition Probabilities 
q      = zeros(2,2);
q(1,1) = 1-1/eta;
q(2,1) = gamma*(1-q(1,1))/(1-gamma);
q(1,2) = 1-q(1,1);
q(2,2) = 1-q(2,1);

% Compute Expected Lifetime Income
e = (eye(2,2)-delta*q)\y;

% Compute Stationary Disbritution of Employment Expected Employment State Durations
p   = markov(q);
eta = 1./(1-diag(q));

% Print State Transition Probabilities
fprintf('\n\n')
fprintf('State Transition Probabilities\n')
fprintf('           Unemployed   Employed\n')
fprintf('Unemployed    %7.3f  %9.3f\n',q(1,:)')
fprintf('Employed      %7.3f  %9.3f\n',q(2,:)')

% Print Stationary Disbritution and Expected State Durations
fprintf('\n\n')
fprintf('Stationary Distribution and Expected State Durations\n')
fprintf('           Unemployed   Employed\n')
fprintf('Proability    %7.3f  %9.3f\n',p(:))
fprintf('Duration      %7.3f  %9.3f\n',eta(:))


%% SIMULATION

% Enter number of replications
n = 100000;

% Perform simulations
iinit = 1;
isim = [iinit; zeros(n-1,1)];
for j=2:n
  isim(j) = markovsim(isim(j-1),q);
end

% Compute visit probabilities
p = [sum(isim==1)/n sum(isim==2)/n];

% Compute average length of runs
run = zeros(n,2);
run(1,iinit) = 1;
for j=2:n
  if isim(j)==1
    if isim(j-1)== 1
      run(j,1) = run(j-1,1)+1;
    end
  end
  if isim(j)==2
    if isim(j-1)== 2
      run(j,2) = run(j-1,2)+1;
    end
  end
end
rcount = sum(~run==0);
eta = sum(run)./rcount;

% Print Simulated Stationary Distribution and Expected State Durations
fprintf('\n\n')
fprintf('Simulated Stationary Distribution and Expected State Durations\n')
fprintf('           Unemployed   Employed\n')
fprintf('Proability    %7.3f  %9.3f\n',p(:))
fprintf('Duration      %7.3f  %9.3f\n',eta(:))