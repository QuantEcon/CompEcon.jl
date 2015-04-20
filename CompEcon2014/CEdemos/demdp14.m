function demdp14alt

%% DEMDP14alt Livestock Feeding Model - Solution of Bellman Equation
%
% Farmer must decide how much to feed his livestock over a finite horizon.
% 
%  States
%      s       weight of livestock at beginning of period
%  Actions
%      x       weight gain current period
%  Parameters
%      alpha   weight gain function parameter
%      beta    weight gain function parameter
%      k       unit cost of feed
%      p       price of livestock per pound
%      N       number of feeding periods
%      s1      initial livestock weight
%      delta   discount factor

% Preliminary tasks
demosetup(mfilename)


%% FORMULATION

% Model Parameters
alpha = 0.1;                           % weight gain function parameter
beta  = 2.0;                           % weight gain function parameter
k     = 0.4;                           % unit cost of feed
p     = 1.0;                           % price of livestock per pound
T     = 6;                             % number of feeding periods
s1    = 0.4;                           % initial livestock weight
delta = 0.95;                          % discount factor


%% SOLUTION - SOLVE EULER EQUATION AS ROOTFINDING PROBLEM

% Solve Euler Conditions
tic
X = [s1*ones(T+1,1);p*ones(T+1,1)];
X = broyden(@f,X,alpha,beta,k,s1,p,delta);
toc
s = X(1:T+1);
l = X(T+2:end);

% Print Terminal Weight
fprintf('Terminal Weight = %5.2f\n',s(end))

% Plot State Path
figure
plot(1:T+1,s)
title('State Path')
xlabel('Period')
ylabel('Livestock Weight')
xlim([1 T+1])

% Plot Action Path
figure
plot(1:T,diff(s))
title('State Path')
xlabel('Period')
ylabel('Weight Gain')
xlim([1 T])

% Plot Shadow Price
figure
plot(1:T+1,l)
title('Shadow Price')
xlabel('Period')
ylabel('Price')
xlim([1 T+1])


%% SOLUTION - SOLVE BELLMAN EQUATION 

% Model Structure
model.horizon = T;                     % number of decision periods
model.func = @func;                    % model functions
model.discount = delta;                % discount factor
model.params = {alpha beta k};         % other parameters
model.ds = 1;                          % dimension of continuous state
model.dx = 1;                          % dimension of continuous action
model.ni = 0;                          % number of discrete states
model.nj = 0;                          % number of discrete actions

% Approximation Structure
n    = 50;                             % number of collocation nodes
smax = 4.0;                            % maximum state
basis = fundefn('spli',n,s1,smax);     % basis functions
snodes = funnode(basis);               % collocaton nodes

% Intialize Policy and Set Post-Terminal Value Function
vterm = p*snodes;

% Solve Bellman Equations
[c,s,v,x] = dpsolve(model,basis,vterm);

% Plot Optimal Policy
figure
plot(s,x)
title('Optimal Weight Gain Policy')
legend('t=1','t=2','t=3','t=4','t=5','t=6')
xlabel('Livestock Weight')
ylabel('Weight Gain')

% Plot Value Function
figure
plot(s,v)
legend('t=1','t=2','t=3','t=4','t=5','t=6','t=7','Location','NW')
title('Value Function')
xlabel('Livestock Weight')
ylabel('Value')

% Plot Shadow Price Function
figure
p = funeval(c,basis,s,1);
plot(s,p)
legend('t=1','t=2','t=3','t=4','t=5','t=6','t=7')
title('Shadow Price Function')
xlabel('Livestock Weight')
ylabel('Price')


%% SIMULATION

% Simulate Model
[ssim,xsim] = dpsimul(model,basis,T,s1,1,s,v,x);

% % Compute and Print Terminal Weight
% fprintf('Terminal Weight = %5.2f\n',ssim(end))
% 
% % Plot State Path
% figure
% plot(1:T+1,ssim)
% title('State Path')
% xlabel('Period')
% ylabel('Livestock Weight')
% xlim([1 T+1])
% 
% % Plot Policy Path
% figure
% plot(1:T,xsim)
% title('Policy Path')
% xlabel('Period')
% ylabel('Weight Gain')
% xlim([1 T])


%% Save Plots as EPS Files
printfigures(mfilename,6)

  
%% Function File
function y = f(X,alpha,beta,k,s0,p,delta)
N = length(X)/2;
s = X(1:N);
l = X(N+1:2*N);
y = [s(1)-s0; ...
    -k*beta*(s(2:N) - s(1:N-1) + alpha*s(1:N-1)).^(beta-1) + delta*l(2:N); ...
    -k*alpha*beta*(s(2:N) - s(1:N-1) + alpha*s(1:N-1)).^(beta-1) + delta*l(2:N) - l(1:N-1); ...
    l(N) - p];

  
%% DPSOLVE FUNCTION FILE
%
%    User-supplied function called by dpsolve that returns the bound,
%    reward, and continuous state transition function values and
%    derivatives with respect to the continuous action x at an arbitrary
%    number ns of states and actions aaccording to the format
%      [out1,out2,out3] = func(flag,s,x,i,j,in,e,<params>)
%    For further information regarding input format, type > help dpsolve or
%    reference the DPSOLVE Handbook.

function [out1,out2,out3] = func(flag,s,x,i,j,in,e,alpha,beta,k)

ns = length(s);
switch flag
case 'b'
  out1 = zeros(ns,1);
  out2 = inf*ones(ns,1);
  out3 = [];
case 'f'
  out1 = -k*(x+alpha*s).^beta;
  out2 = -k*beta*(x+alpha*s).^(beta-1);
  out3 = -k*beta*(beta-1)*(x+alpha*s).^(beta-2);
case 'g'
  out1 = s + x;
  out2 = ones(ns,1);
  out3 = zeros(ns,1);
end  
%%