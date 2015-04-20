function demdp13

%% DEMDP13 Inventory Management Model
%
% Profit maximizing entrepeneur must decide how much to produce and how 
% much inventory to hold.
%
% States
%     s1      market price
%     s2      stock of inventory
% Actions
%     x1      quantity produced
%     x2      inventory carryout
% Parameters
%     c       production cost function parameters
%     k       inventory holding cost function parameters
%     pbar	  long-run mean price
%     rho     mean-reversion coefficient
%     sigma   standard deviation of price shocks
%     delta   discount factor

% Preliminary tasks
demosetup(mfilename)


%% FORMULATION

% Model Parameters
c     = [0.5 0.1];                      % production cost function parameters
k     = [0.1 0.1];                      % inventory holding cost function parameters
pbar  = 1.0;                            % long-run mean price
rho   = 0.5;                            % mean-reversion coefficient
sigma = 0.2;                            % standard deviation of price shocks
delta = 0.9;                            % discount factor

% Continuous State Shock Distribution
m   = 3;                                % number of price shocks
[e,w] = qnwnorm(m,0,sigma^2);           % price shocks and probabilities

% Model Structure
model.func = @func;                     % model functions
model.params = {c k pbar rho};          % function parameters
model.discount = delta;                 % discount factor
model.ds = 2;                           % dimension of continuous state
model.dx = 2;                           % dimension of continuous action
model.ni = 0;                           % number of discrete states
model.nj = 0;                           % number of discrete actions
model.e  = e;                           % continuous state shocks
model.w  = w;                           % continuous state shock probabilities

% Approximation Structure
n = [15 60];                           	% number of collocation node coordinates, per dimension
smin = [0.2 0.0];                     	% minimum states
smax = [2.5 3.0];                     	% maximum states
basis = fundefn('spli',n,smin,smax);    % basis functions

% Deterministic Steady-State
sstar = [pbar 0];                       % deterministic steady-state state
xstar = [(pbar-c(1))/c(2) 0];           % deterministic steady-state action

% Check Model Derivatives
dpcheck(model,sstar,xstar);


%% SOLUTION

% Solve Bellman Equation
[c,s,v,x,resid] = dpsolve(model,basis);

% Reshape Output for Plotting
nr = n*10;
s1 = reshape(s(:,1),nr);
s2 = reshape(s(:,2),nr);
v = reshape(v,nr);
x1 = reshape(x(:,1),nr);
x2 = reshape(x(:,2),nr);
resid = reshape(resid,nr);

% Compute Shadow Prices
p1 = funeval(c,basis,s,[1 0]);
p1 = reshape(p1,nr);

% Plot Optimal Policy 1
figure
surf(s1,s2,x1,'FaceColor','interp','EdgeColor','interp')
title('Optimal Production Policy')
xlabel('Market Price')
ylabel('Inventory')
zlabel('Production')

% Plot Optimal Policy 2
figure
surf(s1,s2,x2,'FaceColor','interp','EdgeColor','interp')
title('Optimal Inventory Policy')
xlabel('Market Price')
ylabel('Inventory')
zlabel('Carryout')

% Plot Value Function
figure
surf(s1,s2,v,'FaceColor','interp','EdgeColor','interp')
title('Value Function')
xlabel('Market Price')
ylabel('Inventory')
zlabel('Value')

% Plot Shadow Price Function 1
figure
surf(s1,s2,p1,'FaceColor','interp','EdgeColor','interp')
title('Shadow Price of Inventories')
xlabel('Market Price')
ylabel('Inventory')
zlabel('Price')

% Plot Residual
figure
surf(s1,s2,resid,'FaceColor','interp','EdgeColor','interp')
title('Bellman Equation Residual')
xlabel('Market Price')
ylabel('Inventory')
zlabel('Residual')


%% SIMULATION

% Simulate Model
rand('seed',0.945);
nper = 26;
nrep = 50000;
sinit = [pbar*ones(nrep,1) zeros(nrep,1)];
iinit = ones(nrep,1);
[ssim,xsim] = dpsimul(model,basis,nper,sinit,iinit,s,v,x);
s1sim = ssim(:,:,1);
s2sim = ssim(:,:,2);
x1sim = xsim(:,:,1);
x2sim = xsim(:,:,2);

% Compute Ergodic Moments
s1avg = mean(s1sim(:)); 
s2avg = mean(s2sim(:)); 
x1avg = mean(x1sim(:)); 
s1std = std(s1sim(:)); 
s2std = std(s2sim(:)); 
x1std = std(x1sim(:)); 

% Print Ergodic Moments
fprintf('Ergodic Means\n') 
fprintf('   Market Price = %5.2f\n'    ,s1avg)
fprintf('   Inventory    = %5.2f\n'    ,s2avg)
fprintf('   Production   = %5.2f\n\n'  ,x1avg)
fprintf('Ergodic Standard Deviations\n') 
fprintf('   Market Price = %5.2f\n'    ,s1std)
fprintf('   Inventory    = %5.2f\n'    ,s2std)
fprintf('   Production   = %5.2f\n\n'  ,x1std)

% Plot Simulated and Expected State Paths 1
figure
hold on
plot(0:nper-1,s1sim(1:3,:))
plot(0:nper-1,mean(s1sim),'k')
plot(nper-1,s1avg,'k*')
title('Simulated and Expected State Paths')
xlabel('Period')
ylabel('Market Price')

% Plot Simulated and Expected State Paths 2
figure
hold on
plot(0:nper-1,s2sim(1:3,:))
plot(0:nper-1,mean(s2sim),'k')
plot(nper-1,s2avg,'k*')
title('Simulated and Expected State Paths')
xlabel('Period')
ylabel('Inventory')

% Plot Simulated and Expected Policy Paths
figure
hold on
plot(0:nper-1,squeeze(x1sim(1:3,:)))
plot(0:nper-1,mean(x1sim),'k')
plot(nper-1,x1avg,'k*')
title('Simulated and Expected Policy Paths')
xlabel('Period')
ylabel('Production')


%% Save Plots as EPS Files
printfigures(mfilename,8)


%% DPSOLVE FUNCTION FILE
%
%    User-supplied function called by dpsolve that returns the bound,
%    reward, and continuous state transition function values and
%    derivatives with respect to the continuous action x at an arbitrary
%    number ns of states and actions aaccording to the format
%      [out1,out2,out3] = func(flag,s,x,i,j,in,e,<params>)
%    For further information regarding input format, type > help dpsolve or
%    reference the DPSOLVE Handbook.

function [out1,out2,out3] = func(flag,s,x,i,j,in,e,c,k,pbar,rho)

n = size(s,1);
ds = 2;
dx = 2;

switch flag
case 'b'
  out1 = zeros(size(s));
  out2 = inf*ones(size(s));
  out3 = [];
case 'f'
  out2 = zeros(n,dx);
  out3 = zeros(n,dx,dx);
  out1 = s(:,1).*(s(:,2)+x(:,1)-x(:,2)) ...
        - (c(1)+0.5*c(2)*x(:,1)).*x(:,1) ...
        - (k(1)+0.5*k(2)*x(:,2)).*x(:,2);
  out2(:,1) =  s(:,1) - (c(1)+c(2)*x(:,1));
  out2(:,2) = -s(:,1) - (k(1)+k(2)*x(:,2));
  out3(:,1,1) = -c(2)*ones(n,1);
  out3(:,2,2) = -k(2)*ones(n,1);
case 'g'
  out2 = zeros(n,ds,dx);
  out3 = zeros(n,ds,dx,dx);
  out1 = [pbar+rho*(s(:,1)-pbar)+e x(:,2)];
  out2(:,2,2) = ones(n,1);
end
%%