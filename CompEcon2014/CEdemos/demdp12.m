function demdp12

%% DEMDP12 Production Management Model
%
% Profit maximizing entrepeneur must decide how much to produce, subject to 
% production adjustment costs.
%
% States
%     i       market price (discrete)
%     s       lagged production (continuous)
% Actions
%     x       current production
% Parameters
%     alpha   marginal adjustment cost
%     beta    marginal production cost parameters
%     pbar    long-run average market price
%     sigma   market price shock standard deviation
%     delta   discount factor

% Preliminary tasks
demosetup(mfilename)


%% FORMULATION

% Model Parameters
alpha = 0.01;                               % marginal adjustment cost
beta  = [0.8 0.03];                         % marginal production cost parameters
pbar  = 1.0;                                % long-run average market price
sigma = 0.2;                                % market price shock standard deviation
delta = 0.9;                                % discount factor

% Continuous State Shock Distribution
m = 3;                                      % number of market price shocks
mu = log(pbar)-sigma^2/2;                   % mean log price
[p,w] = qnwlogn(m,mu,sigma^2);              % market price shocks and probabilities
q = w(:,ones(1,m))';                        % market price transition probabilities

% Model Structure
model.func = @func;                         % model functions
model.params = {alpha beta p};              % function parameters
model.discount = delta;                     % discount factor
model.ds = 1;                               % dimension of continuous state
model.dx = 1;                               % dimension of continuous action
model.ni = m;                               % number of discrete states
model.nj = 0;                               % number of discrete actions
model.q  = q;                               % discrete state transition probabilities

% Approximation Structure
n = 50;                                     % number of collocation nodes
smin =  0;                                  % minimum state
smax = 20;                                  % maximum state
basis = fundefn('spli',n,smin,smax);        % basis functions

% Deterministic Steady-State
sstar = (pbar-beta(1))/beta(2);             % deterministic steady-state state

% Check Model Derivatives
dpcheck(model,sstar,sstar);


%% SOLUTION

% Solve Bellman Equation
[c,s,v,x,resid] = dpsolve(model,basis);

% Plot Optimal Policy
figure
plot(s,x)
legend('Low Price','Average Price','High Price')
legend('Location','Best')
legend('boxoff')
title('Optimal Production Policy')
xlabel('Lagged Production')
ylabel('Production')

% Plot Value Function
figure
plot(s,v)
legend('Low Price','Average Price','High Price')
legend('Location','Best')
legend('boxoff')
title('Value Function')
xlabel('Lagged Production')
ylabel('Value')

% Plot Shadow Price Function
figure
lambda = funeval(c,basis,s,1);
plot(s,lambda)
legend('Low Price','Average Price','High Price')
legend('Location','Best')
legend('boxoff')
title('Shadow Price of Lagged Production')
xlabel('Lagged Production')
ylabel('Shadow Price')

% Plot Residual
figure
hold on
plot(s,resid)
legend('Low Price','Average Price','High Price')
legend('Location','Best')
legend('boxoff')
plot(s,0*resid,'k--','LineWidth',2)
title('Bellman Equation Residual')
xlabel('Lagged Production')
ylabel('Residual')


%% SIMULATION

% Simulate Model
rand('seed',0.945);
nper = 26;
nrep = 10000;
sinit = smin(ones(nrep,1),:);
iinit = 2*ones(nrep,1);
[ssim,xsim,isim] = dpsimul(model,basis,nper,sinit,iinit,s,v,x);
psim = p(isim);

% Compute Ergodic Moments
savg = mean(ssim(:)); 
pavg = mean(psim(:)); 
sstd = std(ssim(:)); 
pstd = std(psim(:)); 

% Print Ergodic Mean and Standard Deviation
fprintf('Deterministic Steady-State\n') 
fprintf('   Price       = %5.2f\n'  ,pbar)
fprintf('   Production  = %5.2f\n\n',sstar)
fprintf('Ergodic Means\n') 
fprintf('   Price       = %5.2f\n',  pavg)
fprintf('   Production  = %5.2f\n\n',savg)
fprintf('Ergodic Standard Deviations\n') 
fprintf('   Price       = %5.2f\n'  ,pstd)
fprintf('   Production  = %5.2f\n\n',sstd)

% Plot Simulated and Expected Policy Path
figure
hold on
plot(0:nper-1,xsim(1:3,:))
plot(0:nper-1,mean(xsim),'k')
title('Simulated and Expected Production')
xlabel('Period')
ylabel('Production')


%% Save Plots as EPS Files
printfigures(mfilename,5)


%% DPSOLVE FUNCTION FILE
%
%    User-supplied function called by dpsolve that returns the bound,
%    reward, and continuous state transition function values and
%    derivatives with respect to the continuous action x at an arbitrary
%    number ns of states and actions aaccording to the format
%      [out1,out2,out3] = func(flag,s,x,i,j,in,e,<params>)
%    For further information regarding input format, type > help dpsolve or
%    reference the DPSOLVE Handbook.

function [out1,out2,out3] = func(flag,s,q,i,j,in,e,alpha,beta,p)

n = length(s);
p = p(i);
l = s;
switch flag
case 'b'
  out1 = zeros(n,1);
  out2 = inf*ones(n,1);
  out3 = [];
case 'f'
  out1 = p.*q - (beta(1)*q+0.5*beta(2)*q.^2) - 0.5*alpha*((q-l).^2);
  out2 = p - beta(1) - beta(2)*q - alpha*(q-l);
  out3 = (-beta(2)-alpha)*ones(n,1);
case 'g'
  out1 = q;
  out2 = ones(n,1);
  out3 = zeros(n,1);
end
%%