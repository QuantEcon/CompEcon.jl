function demdp09

%% DEMDP09 Private Non-Renewable Resource Model
%
%  Profit maximizing mine owner must decide how much ore to extract.
%
% States
%     s       ore stock
% Actions
%     q       ore extracted and sold
% Parameters
%     a       demand function parameters
%     b       cost function parameters
%     delta   discount factor

% Preliminary tasks
demosetup(mfilename)


%% FORMULATION
 
% Model Parameters
a = [5 0.8];                            % demand function parameters
b = [7 1.0];                            % cost function parameters
delta = 0.9;                            % discount factor

% Model Structureoeeeeadfasd
model.func = @func;                     % model functions
model.params = {a b};                   % function parameters
model.discount = delta;                 % discount factor
model.ds = 1;                           % dimension of continuous state
model.dx = 1;                           % dimension of continuous action
model.ni = 0;                           % number of discrete states
model.nj = 0;                           % number of discrete actions

% Approximation Structure
n    = 101;                             % number of collocation nodes
smin =   0;                             % minimum state
smax =  10;                             % maximum state
basis = fundefn('spli',n,smin,smax);    % basis functions

% Check Model Derivatives
dpcheck(model,smax,0);


%% SOLUTION

% Solve Bellman Equation
[c,s,v,q,resid] = dpsolve(model,basis);

% Compute and print abandonment point
sstar = (b(1)-a(1))/b(2);
fprintf('Abandonment Point = %5.2f\n',sstar)

% Plot Optimal Policy
figure
hold on
plot(s,q)
title('Optimal Extraction')
xlabel('Ore Stock')
ylabel('Ore Extracted')

% Plot Value Function
figure
plot(s,v)
title('Value Function')
xlabel('Ore Stock')
ylabel('Value')

% Plot Shadow Price Function
figure
p = funeval(c,basis,s,1);
plot(s,p)
title('Shadow Price Function')
xlabel('Ore Stock')
ylabel('Shadow Price')

% Plot Residual
figure
hold on
plot(s,resid)
plot(s,0*s,'k--')
title('Bellman Equation Residual')
xlabel('Ore Stock')
ylabel('Residual')


%% SIMULATION

% Simulate Model
nper = 21; 
sinit = smax;
iinit = 1;
[ssim,qsim] = dpsimul(model,basis,nper,sinit,iinit,s,v,q);

% Plot State and Policy Paths
figure
hold on
plot(0:nper-1,ssim,0:nper-1,qsim)
title('State and Policy Paths')
legend('Stock','Extraction')
legend('Location','Best')
legend boxoff
plot([0 nper-1],[sstar sstar],'k--')
xlabel('Period')
ylabel('Stock/Extraction')


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

function [out1,out2,out3] = func(flag,s,q,i,j,in,e,a,b)

switch flag
case 'b'
  out1 = zeros(size(s));
  out2 = s;
  out3 = [];
case 'f'
  out1 = (a(1)-b(1)+b(2)*s).*q - (a(2)+b(2)/2).*q.^2;
  out2 = (a(1)-b(1)+b(2)*s) - 2*(a(2)+b(2)/2).*q;
  out3 = -2*(a(2)+b(2)/2)*ones(size(s));
case 'g'
  out1 = s-q;
  out2 = -ones(size(s));
  out3 = zeros(size(s));
end
%%