function demdp08

%% DEMDP08 Public Renewable Resource Model
%
% Welfare maximizing social planner must decide how much of a renewable 
% resource to harvest.
%
% States
%     s       quantity of stock available
% Actions
%     q       quantity of stock harvested
% Parameters
%     alpha   growth function parameter
%     beta    growth function parameter
%     gamma   relative risk aversion
%     kappa   unit cost of harvest
%     delta   discount factor

% Preliminary tasks
demosetup(mfilename)


%% FORMULATION
  
% Model Parameters
alpha = 4.0;                                    % growth function parameter
beta  = 1.0;                                  	% growth function parameter
gamma = 0.5;                                   	% relative risk aversion
kappa = 0.2;                                    % unit cost of harvest
delta = 0.9;                                  	% discount factor
    
% Model Structure
model.func = @func;                             % model functions
model.params = {alpha beta gamma kappa};        % function parameters
model.discount = delta;                         % discount factor
model.ds = 1;                                   % dimension of continuous state
model.dx = 1;                                   % dimension of continuous action
model.ni = 0;                                   % number of discrete states
model.nj = 0;                                   % number of discrete actions

% Approximation Structure
n    = 8;                                     	% number of collocation nodes
smin = 6;                                      	% minimum state
smax = 9;                                       % maximum state
basis = fundefn('cheb',n,smin,smax);            % basis functions
  
% Steady-State
sstar = (alpha^2-1/delta^2)/(2*beta);           % steady-state stock
qstar = sstar - (delta*alpha-1)/(delta*beta); 	% steady-state action

% Print Steady-States
fprintf('Steady States\n') 
fprintf('   Stock   = %5.2f\n'  ,sstar)
fprintf('   Harvest = %5.2f\n'  ,qstar)

% Check Model Derivatives
dpcheck(model,sstar,qstar);


%% SOLUTION

% Solve Bellman Equation
[c,s,v,q,resid] = dpsolve(model,basis);

% Plot Optimal Policy
figure
plot(s,q)
title('Optimal Harvest Policy')
xlabel('Stock')
ylabel('Harvest')

% Plot Value Function
figure
plot(s,v)
title('Value Function')
xlabel('Stock')
ylabel('Value')

% Plot Shadow Price Function
figure
p = funeval(c,basis,s,1);
plot(s,p)
title('Shadow Price Function')
xlabel('Stock')
ylabel('Shadow Price')

% Plot Residual
figure
hold on
plot(s,resid)
plot(s,0*s,'k--')
title('Bellman Equation Residual')
xlabel('Stock')
ylabel('Residual')


%% SIMULATION

% Simulate Model
nper = 16; 
sinit = smin;
iinit = 1;
[ssim,qsim] = dpsimul(model,basis,nper,sinit,iinit,s,v,q);

% Plot State and Policy Paths
figure
hold on
plot(0:nper-1,ssim,0:nper-1,qsim)
title('State and Policy Paths')
legend('Stock','Harvest')
legend('Location','Best')
legend boxoff
plot(nper-1,sstar,'*',nper-1,qstar,'*')
xlabel('Period')
ylabel('Stock/Harvest')


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

function [out1,out2,out3] = func(flag,s,q,i,j,in,e,alpha,beta,gamma,kappa)

switch flag
case 'b'
  out1 = zeros(size(s));                  
  out2 = s;        
  out3 = [];                      
case 'f'
  out1 = (q.^(1-gamma))/(1-gamma)-kappa*q; 
  out2 = q.^(-gamma)-kappa;                
  out3 = -gamma*q.^(-gamma-1);            
case 'g'
  out1 = alpha*(s-q) - 0.5*beta*(s-q).^2; 
  out2 = -alpha + beta*(s-q);             
  out3 = -beta*ones(size(s));  
end
%%