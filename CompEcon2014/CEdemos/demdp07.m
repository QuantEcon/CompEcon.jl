function demdp07

%% DEMDP07 Stochastic Optimal Economic Growth Model
%
% Welfare maximizing social planner must decide how much society should
% consume and invest. Unlike the deterministic model, this model allows
% arbitrary constant relative risk aversion, capital depreciation, and 
% stochastic production shocks.  It lacks a known closed-form solution.
%
% States
%     s       stock of wealth
% Actions
%     k       capital investment
% Parameters
%     alpha   relative risk aversion
%     beta    capital production elasticity
%     gamma   capital survival rate
%     sigma   production shock volatility
%     delta   discount factor

% Preliminary tasks
demosetup(mfilename)


%% FORMULATION
  
% Model Parameters
alpha = 0.2;                            % relative risk aversion
beta  = 0.5;                            % capital production elasticity
gamma = 0.9;                          	% capital survival rate
sigma = 0.1;                            % production shock volatility
delta = 0.9;                            % discount factor

% Continuous State Shock Distribution
m = 5;                                 	% number of production shocks
[e,w] = qnwlogn(m,-sigma^2/2,sigma^2); 	% production shocks and probabilities

% Model Structure
model.func = @func;                     % model functions
model.params = {alpha beta gamma};	    % function parameters
model.discount = delta;                 % discount factor
model.ds = 1;                           % dimension of continuous state
model.dx = 1;                           % dimension of continuous action
model.ni = 0;                           % number of discrete states
model.nj = 0;                           % number of discrete actions
model.e  = e;                           % continuous state shocks
model.w  = w;                           % continuous state shock probabilities

% Approximation Structure
n     = 10;                             % number of collocation nodes
smin  =  5;                             % minimum wealth
smax  = 10;                             % maximum wealth
basis = fundefn('cheb',n,smin,smax);    % basis functions
  
% Deterministic Steady-State
kstar = ((1-delta*gamma)/(delta*beta))^(1/(beta-1)); % determistic steady-state capital investment
sstar = gamma*kstar + kstar^beta;      	% deterministic steady-state wealth

% Check Model Derivatives
dpcheck(model,sstar,kstar);


%% SOLUTION

% Solve Bellman Equation
[c,s,v,k,resid] = dpsolve(model,basis);
 
% Plot Optimal Policy
figure
plot(s,k)
title('Optimal Investment Policy')
xlabel('Wealth')
ylabel('Investment')

% Plot Value Function
figure
plot(s,v)
title('Value Function')
xlabel('Wealth')
ylabel('Value')

% Plot Shadow Price Function
figure
pr = funeval(c,basis,s,1);
plot(s,pr)
title('Shadow Price Function')
xlabel('Wealth')
ylabel('Shadow Price')

% Plot Residual
figure
plot(s,resid,s,0*s,'k--')
title('Bellman Equation Residual')
xlabel('Wealth')
ylabel('Residual')


%% SIMULATION

% Simulate Model
rng('default')
nper = 21; 
nrep = 50000;
sinit = smin*ones(nrep,1);
iinit = ones(nrep,1);
[ssim,ksim] = dpsimul(model,basis,nper,sinit,iinit,s,v,k);

% Plot Simulated State Path
figure
hold on
plot(0:nper-1,ssim(1:3,:))
plot(0:nper-1,mean(ssim),'k')
title('Simulated and Expected Wealth')
xlabel('Period')
ylabel('Wealth')

% Plot Simulated Policy Path
figure
hold on
plot(0:nper-1,ksim(1:3,:))
plot(0:nper-1,mean(ksim),'k')
title('Simulated and Expected Investment')
xlabel('Period')
ylabel('Investment')

% Compute Ergodic Moments
ssim = ssim(:,nper);
ksim = ksim(:,nper);
savg = mean(ssim); 
kavg = mean(ksim); 
sstd = std(ssim); 
kstd = std(ksim); 

% Print Steady-State and Ergodic Moments
fprintf('Deterministic Steady-State\n') 
fprintf('   Wealth       = %5.2f\n'    ,sstar)
fprintf('   Investment   = %5.2f\n\n'  ,kstar)
fprintf('Ergodic Means\n') 
fprintf('   Wealth       = %5.2f\n'    ,savg)
fprintf('   Investment   = %5.2f\n\n'  ,kavg)
fprintf('Ergodic Standard Deviations\n') 
fprintf('   Wealth       = %5.2f\n'    ,sstd)
fprintf('   Investment   = %5.2f\n\n'  ,kstd)

% Compute and Plot Ergodic Wealth Distribution
[qq,ss] = ksdensity(ssim(:),'support','positive');
figure
plot(ss,qq)
title('Ergodic Wealth Distribution')
xlabel('Wealth')
ylabel('Probability')
xlim([5 10])  


%% Save Plots as EPS Files
printfigures(mfilename,7)


%% DPSOLVE FUNCTION FILE
%
%    User-supplied function called by dpsolve that returns the bound,
%    reward, and continuous state transition function values and
%    derivatives with respect to the continuous action x at an arbitrary
%    number ns of states and actions aaccording to the format
%      [out1,out2,out3] = func(flag,s,x,i,j,in,e,<params>)
%    For further information regarding input format, type > help dpsolve or
%    reference the DPSOLVE Handbook.

function [out1,out2,out3] = func(flag,s,k,i,j,in,e,alpha,beta,gamma)

switch flag
case 'b'
  out1 = zeros(size(s));  
  out2 = 0.99*s;  
  out3 = [];
case 'f'
  out1 = ((s-k).^(1-alpha))/(1-alpha); 
  out2 = -(s-k).^(-alpha);  
  out3 = -alpha*(s-k).^(-alpha-1);
case 'g'
  out1 = gamma*k + e.*k.^beta;          
  out2 = gamma + beta*e.*k.^(beta-1); 
  out3 = (beta-1)*beta*e.*k.^(beta-2);
end 
%%