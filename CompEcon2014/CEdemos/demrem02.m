function demrem02

%% DEMREM02 Commodity Storage Model
%
% Solves standard rational expectations commodity storage model.
%
% States
%     s         supply
% Actions
%     x         carryout
% Parameters
%     gamma     inverse demand elasticity
%     xbar      storage capacity
%     k         unit storage cost
%     sigma     production volatility
%     delta     discount factor

% Preliminary tasks
demosetup(mfilename)


%% FORMULATION

% Model Parameters
gamma = 0.5;                              % inverse demand elasticity
xbar  = 1.0;                              % storage capacity
k     = 0.05;                             % unit storage cost
sigma = 0.3;                              % production  volatility
delta = 0.95;                             % discount factor

% Shock Distribution
m = 5;                                    % number of production shocks
[y,w] = qnwlogn(m,-0.5*sigma^2,sigma^2);  % log normal nodes and weights

% Model Structure
model.func = @func;                       % model functions
model.params = {delta,gamma,xbar,k};      % other parameters
model.ds = 1;                             % dimension of state s
model.dx = 1;                             % dimension of response x
model.e = y;                              % production levels
model.w = w;                              % probabilities

% Approximation Structure
n      = 150;                             % degree of approximation
smin   = 0.3;                             % minimum supply
smax   = 5.0;                             % maximum supply
basis  = fundefn('spli',n,smin,smax);     % basis functions
s = funnode(basis);                       % collocation nodes


%% SOLUTION

% Compute Rational Expectations Equilibrium Directly
tic
x = zeros(n,1);
c = zeros(n,1);
for it=1:100
  cold = c;
  f = -(s-x).^(-gamma)-k;
  d = -gamma*(s-x).^(-gamma-1);
  for j=1:m
    sn = x + y(j);
    f = f + w(j)*delta*funeval(c,basis,sn);
    d = d + w(j)*delta*funeval(c,basis,sn,1);
  end
  x = x + min(max(-f./d,-x),xbar-x);
  p = (s-x).^(-gamma);
  c = funfitxy(basis,s,p);
  change = norm(c-cold,inf);
  fprintf ('%4i %10.1e\n',it,change)
  if change<1.e-10, break, end
end
toc

% Compute Rational Expectations Equilibrium Using remsolve
[c,s,x,f,resid] = remsolve(model,basis);

% Plot Equilibrium Response 1
figure
plot(s,x)
title('Equilibrium Carryout')
xlabel('Supply')
ylabel('Carryout')
xlim([smin smax])
ylim([0 1.0])

% Plot Equilibrium Response 2
figure
hold on
plot(s,(s-x).^-gamma,'r')
plot(s,s.^-gamma)
legend('Total Demand','Consumption Demand')
legend boxoff
% plot([0 smax],[pstar pstar],'k:')
title('Equilibrium Price')
xlabel('Quantity')
ylabel('Price')

% Plot Expected Arbitrage Benefit
figure
hold on
plot(s,f)
plot(s,0*f,'k--','LineWidth',2)
title('Expected Arbitrage Profit')
xlabel('Supply')
ylabel('Profit')
xlim([smin smax])

% Plot Response Residual
figure
hold on
plot(s,resid)
plot(s,0*resid,'k--','LineWidth',2)
title('Response Residual')
xlabel('Supply')
ylabel('Residual')
xlim([smin smax])


%% SIMULATION

% Intializations
nper  = 20;
nrep  = 10000;
sinit = ones(nrep,1);
ssim  = zeros(nrep,nper+1);
xsim  = zeros(nrep,nper+1);

% Simulate Model Not Using remsimul
tic
rand('seed',0.945);
ss = sinit;
for ip=1:nper+1
  xx = interp1(s,x,ss,'cubic');
  ssim(:,ip,:) = ss;
  xsim(:,ip,:) = xx;
  if ip<nper+1
    ss = xx + y(discrand(nrep,w),:);
  end
end
psim = (ssim-xsim).^(-gamma);
toc

% Simulate Model Using remsimul
tic
rand('seed',0.945);
[ssim,xsim] = remsimul(model,basis,nper,sinit,s,x);
psim = (ssim-xsim).^-gamma;
toc

% Plot Simulated and Expected State
figure
plot(0:nper,ssim(1:3,:),0:nper,mean(ssim),'k')
title('Simulated and Expected Supply')
xlabel('Period')
ylabel('Supply')

% Plot Simulated and Expected Response 1
figure
plot(0:nper,xsim(1:3,:),0:nper,mean(xsim),'k')
title('Simulated and Expected Carryout')
xlabel('Period')
ylabel('Carryout')

% Plot Simulated and Expected Response 2
figure
plot(0:nper,psim(1:3,:),0:nper,mean(psim),'k')
title('Simulated and Expected Price')
xlabel('Period')
ylabel('Price')

% Compute Ergodic Moments
savg = mean(ssim(:)); 
xavg = mean(xsim(:)); 
pavg = mean(psim(:)); 
sstd = std(ssim(:)); 
xstd = std(xsim(:)); 
pstd = std(psim(:)); 

% Print Ergodic Moments
fprintf('Ergodic Means\n') 
fprintf('   Supply   = %5.2f\n'    ,savg)
fprintf('   Carryout = %5.2f\n'    ,xavg)
fprintf('   Price    = %5.2f\n\n'  ,pavg)
fprintf('Ergodic Standard Deviations\n') 
fprintf('   Supply   = %5.2f\n'    ,sstd)
fprintf('   Carryout = %5.2f\n'    ,xstd)
fprintf('   Price    = %5.2f\n\n'  ,pstd)

% Compute and Plot Ergodic State Distribution
[qq,ss] = ksdensity(ssim(:),'support','positive','bandwidth',0.15);
figure
plot(ss,qq)
title('Ergodic Supply Distribution')
xlabel('Supply')
ylabel('Probability')
xlim([0 4])  

% Compute and Plot Ergodic Price Distribution
[qq,pp] = ksdensity(psim(:),'support','positive','bandwidth',0.15);
figure
plot(pp,qq)
title('Ergodic Price Distribution')
xlabel('Price')
ylabel('Probability')
xlim([0 3])  


%% Save Plots as EPS Files
printfigures(mfilename,9)


%% Function File
%
%   Called by remsolve
%
%   A user-supplied function that returns the bound, arbitrage, and state 
%   transitions, and and their first derivatives with respect to the 
%   response x, at an arbitrary number ns of states s and responses x
%   according to the format
%     [out1,out2,out3,out4] = func(flag,s,x,sn,xn,e,params)
%   Function File Input
%     flag      : flag indicating function to be evaluated 
%     s         : ns.ds states
%     x         : ns.dx responses
%     sn        : ns.ds states next period
%     xn        : ns.dx responses next period
%     e         : ns.de shocks
%     params    : parameters passed to function file
%   Function File Output
%     if flag = 'b', returns lower and upper bounds on response x
%       out1    : ns.dx lower bounds on response x
%       out2    : ns.dx upper bounds on response x
%       out3    : empty
%       out4    : empty
%     if flag = 'f', returns marginal arbitrage profit and derivatives
%       out1    : ns.dx    marginal arbitrage profit
%       out2    : ns.dx.dx first derivative of f with respect to x
%       out3    : ns.dx.ds first derivative of f with respect to sn
%       out4    : ns.dx.dx first derivative of f with respect to xn
%     if flag = 'g', returns state transition and derivatives
%       out1    : ns.ds    state next period g
%       out2    : ns.ds.dx first derivative of g with respect to x
%       out3    : empty
%       out4    : empty

function [out1,out2,out3,out4] = func(flag,s,x,sn,xn,y,delta,gamma,xbar,k)

ns = length(s);
switch flag
  case 'b'
    out1 = zeros(ns,1);
    out2 = xbar*ones(ns,1);
  case 'f'
    out1 = delta*(sn-xn).^(-gamma)-(s-x).^(-gamma)-k;
    out2 = -gamma*(s-x).^(-gamma-1);
    out3 = -gamma*delta*(sn-xn).^(-gamma-1);
    out4 =  gamma*delta*(sn-xn).^(-gamma-1);
  case 'g'
    out1 = x + y;
    out2 = ones(ns,1);
end
%%