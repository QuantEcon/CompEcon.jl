function demdp06

%% DEMDP06 Deterministic Optimal Economic Growth Model
%
% Welfare maximizing social planner must decide how much society should
% consume and invest.  Model is of special interest because it has a known 
% closed-form solution.
%
% States
%     s       stock of wealth
% Actions
%     k       capital investment
% Parameters
%     beta    capital production elasticity
%     delta   discount factor

% Preliminary tasks
demosetup(mfilename)


%% FORMULATION

% Model Parameters
beta  = 0.5;                            % capital production elasticity
delta = 0.9;                            % discount factor

% Model Structure
model.func = @func;                     % model functions
model.params = {beta};                  % function parameters
model.discount = delta;                 % discount factor
model.ds = 1;                           % dimension of continuous state
model.dx = 1;                           % dimension of continuous action
model.ni = 0;                           % number of discrete states
model.nj = 0;                           % number of discrete actions

% Approximation Structure
n     = 15;                             % number of collocation nodes
smin  = 0.2;                            % minimum wealth
smax  = 1.0;                            % maximum wealth
basis = fundefn('cheb',n,smin,smax);    % basis functions
s = funnode(basis);                     % collocaton nodes

% Steady-State
sstar = (beta*delta)^(beta/(1-beta));   % steady-state wealth
kstar = beta*delta*sstar;               % steady-state capital investment
vstar = log(sstar-kstar)/(1-delta);     % steady-state value
pstar = 1/(sstar*(1-beta*delta));       % steady-state shadow price
b = 1/(1-delta*beta);
  
% Check Model Derivatives
dpcheck(model,sstar,kstar);


%% SOLUTION

% Compute Analytic Solution at Collocation Nodes
vtrue = vstar + b*(log(s)-log(sstar));
ktrue = delta*beta*s;

% Solve Bellman Equation
[c,s,v,k,resid] = dpsolve(model,basis,vtrue,ktrue);
   
% Compute Linear-Quadratic Approximation at Refined State Grid
[vlq,klq,plq] = lqapprox(model,s,sstar,kstar);

% Compute Analytic Solution at Refined State Grid
vtrue = vstar + b*(log(s)-log(sstar));

% Plot Optimal Policy
figure
plot(s,[k klq],[sstar sstar],[-0.2 kstar],'r--')
bullet(sstar,kstar,25,'r')
ftext(sstar,-0.17,'\bf s^*','left','middle',14)
title('Optimal Investment Policy')
legend('Chebychev Collocation','L-Q Approximation','Location','NW')
xlabel('Wealth')
ylabel('Investment')

% Plot Value Function
figure
plot(s,[v vlq],[sstar sstar],[-16 vstar],'r--')
bullet(sstar,vstar,25,'r')
ftext(sstar,-15.8,'\bf s^*','left','middle',14)
title('Value Function')
legend('Chebychev Collocation','L-Q Approximation','Location','NW')
xlabel('Wealth')
ylabel('Value')

% Plot Shadow Price Function
figure
pr = funeval(c,basis,s,1);
plot(s,[pr plq],[sstar sstar],[0 pstar],'r--')
bullet(sstar,pstar,25,'r')
ftext(sstar,0.3,'\bf s^*','left','middle',14)
title('Shadow Price Function')
legend('Chebychev Collocation','L-Q Approximation')
xlabel('Wealth')
ylabel('Shadow Price')

% Plot Chebychev Collocation Residual and Approximation Error
figure
plot(s,[resid v-vtrue],s,0*s,'k--')
title('Chebychev Collocation Residual and Approximation Error')
legend('Residual','Error','Location','SE')
xlabel('Wealth')
ylabel('Residual/Error')

% Plot Linear-Quadratic Approximation Error
figure
plot(s,vlq-vtrue)
title('Linear-Quadratic Approximation Error')
xlabel('Wealth')
ylabel('Error')


%% SIMULATION

% Simulate Model
nper = 21; 
nrep = 1;
sinit = smin*ones(nrep,1);
iinit = ones(nrep,1);
[ssim,ksim] = dpsimul(model,basis,nper,sinit,iinit,s,v,k);

% Plot State and Policy Paths
figure
hold on
plot(0:nper-1,ssim,0:nper-1,ksim)
title('State and Policy Paths')
legend('Wealth','Investment')
legend('Location','Best')
legend boxoff
plot(nper-1,sstar,'*',nper-1,kstar,'*')
xlabel('Period')
ylabel('Wealth/Investment')

% Print Steady-State
fprintf('Steady-State\n') 
fprintf('   Wealth       = %5.2f\n'    ,sstar)
fprintf('   Investment   = %5.2f\n\n'  ,kstar)


%% Save Plots as EPS Files
printfigures(mfilename,6)


%% DPSOLVE FUNCTION FILE
%
%    User-supplied function called by dpsolve that returns the bound,
%    reward, and continuous state transition function values and
%    derivatives with respect to the continuous action x at an arbitrary
%    number ns of states and actions aaccording to the format
%      [out1,out2,out3] = func(flag,s,x,i,j,in,e,<params>)
%    For further information regarding input format, type > help dpsolve or
%    reference the DPSOLVE Handbook.

function [out1,out2,out3] = func(flag,s,k,i,in,j,e,beta)

n = length(s);
switch flag
case 'b'
  out1 = zeros(n,1);  
  out2 = s;  
  out3 = [];
case 'f'
  out1 = log(s-k); 
  out2 = -(s-k).^(-1);
  out3 = -(s-k).^(-2);
case 'g'
  out1 = k.^beta;          
  out2 = beta*k.^(beta-1);
  out3 = (beta-1)*beta*k.^(beta-2); 
end 
%%