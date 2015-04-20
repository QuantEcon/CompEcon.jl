function demdp05

%% DEMDP05 American Put Option Pricing Model
%
% Compute the critical exercise price for an American put option in terms 
% of time to expiration.
%
% States
%     p       underlying asset price
% Actions
%     j       exercize (2) or do not exercize (1) option
% Parameters
%     K       option strike price
%     N       number of periods to expiration
%     mu      mean of log price innovation
%     sigma   standard devitation of log price innovation
%     delta   discount factor

% Preliminary tasks
demosetup(mfilename)


%% FORMULATION

% Continuous Time Discretization
% sigma = 0.2;                            % annual volatility
% T     = 0.5;                            % years to expiration
% K     = 1.0;                            % option strike price
% r     = 0.1;                            % annual interest rate
% N     = 300;                            % number of time intervals
% dt    = T/N;                            % length of time interval
% delta = exp(-r*dt);                     % per period discount factor
% mu    = dt*(r-sigma^2/2);               % mean of log price innovation
  
% Model Parameters
K     = 1.0;                            % option strike price
N     = 300;                            % number of periods to expiration
mu    = 0.0001;                         % mean of log price innovation
sigma = 0.0080;                         % standard devitation of log price innovation
delta = 0.9998;                         % discount factor
  
% Continuous State Shock Distribution
m = 15;                                 % number of price shocks
[e,w] = qnwnorm(m,mu,sigma^2);          % price shocks and probabilities  
  
% Approximation Structure
n    = 500;                             % number of collocation nodes
pmin  = -1;                             % minimum log price
pmax  =  1;                             % maximum log price
basis = fundefn('spli',n,pmin,pmax);    % basis functions  
p     = funnode(basis);                 % collocaton nodes
Phi   = funbase(basis);                 % interpolation matrix


%% SOLUTION
  
% Intialize Value Function
c = zeros(n,1);                         % conditional value function basis coefficients

% Solve Bellman Equation and Compute Critical Exercise Prices
f = @(p,K,delta,c,basis) K-exp(p)-delta*funeval(c,basis,p);
pcrit = broyden(f,0,K,delta,c,basis);
for t=1:N
    v = zeros(n,1);
    for k=1:m
        pnext = p + e(k);
        v = v + w(k)*max(K-exp(pnext),delta*funeval(c,basis,pnext));
    end
    c = Phi\v;
    pcrit = [broyden(f,0,K,delta,c,basis) pcrit];
end

% Print Critical Exercise Price 300 Periods to Expiration
fprintf('Critical Exercise Price 300 Periods to Expiration\n') 
fprintf('   Critical Price  = %5.2f\n\n' ,exp(pcrit(1)))

% Plot Critical Exercise Prices
figure
time = (N:-1:0);
plot(time,exp(pcrit))
title('American Put Option Optimal Exercise Boundary')
xlabel('Periods Remaining Until Expiration')
ylabel('Exercise Price')


%% Save Plots as EPS Files
printfigures(mfilename,1)