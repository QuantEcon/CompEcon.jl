function demrem01

%% DEMREM01 Asset Pricing Model
%
% Solves Lucas-Prescott rational expectations asset pricing model.
%
% States
%     d       asset dividend
% Response
%     p       asset price
% Parameters
%     beta    coefficient of risk aversion
%     dbar    long-run mean dividend
%     gamma   dividend autoregression coefficient
%     sigma   dividend shock standard deviation
%     delta   discount factor

% Preliminary tasks
demosetup(mfilename)


%% FORMULATION

% Model Parameters
beta  = 0.5;                            % coefficient of risk aversion
dbar  = 1.0;                            % long-run mean dividend
gamma = 0.5;                            % dividend autoregression coefficient
sigma = 0.1;                            % dividend shock standard deviation
delta = 0.9;                            % discount factor

% Shock Distribution
m = 5;                                  % number of divident shocks
[e,w] = qnwnorm(m,0,sigma^2);           % normal nodes and weights

% Model Structure
model.func = @func;                     % model functions
model.params = {delta beta dbar gamma}; % function parameters
model.ds = 1;                           % dimension of state s
model.dx = 1;                           % dimension of response x
model.e  = e;                           % state shocks
model.w  = w;                           % state shock probabilities

% Approximation Structure
n      = 25;                            % degree of approximation
dmin   = 0.1;                           % minimum dividend
dmax   = 1.9;                           % maximum dividend
basis  = fundefn('cheb',n,dmin,dmax);   % basis functions
dnode  = funnode(basis);                % collocaton nodes


%% SOLUTION

% Compute Rational Expectations Equilibrium Directly
tic
LHS = diag(dnode.^(-beta))*funbase(basis,dnode);
RHS = 0;
for k=1:m
  dnext = dbar + gamma*(dnode-dbar) + e(k);
  LHS   = LHS - delta*w(k)*diag(dnext.^(-beta))*funbase(basis,dnext);
  RHS   = RHS + delta*w(k)*dnext.^(1-beta);
end
c = LHS\RHS;
p = funeval(c,basis,dnode);

% Compute Residual
dd = nodeunif(10*n,dmin,dmax);
Ef = 0;
for k=1:m
  dnext = dbar + gamma*(dd-dbar) + e(k);
  f     = diag(dnext.^(-beta))*(funeval(c,basis,dnext)+dnext);
  Ef    = Ef + delta*w(k)*f;
end
resid = dd.^(-beta).*funeval(c,basis,dd)-Ef;
toc

% Compute Rational Expectations Equilibrium Using remsolve
[c,d,p,f,resid] = remsolve(model,basis,p);

% Plot Equilibrium Response
figure
plot(d,p)
title('Equilibrium Asset Price')
xlabel('Dividend')
ylabel('Price')
xlim([0.1 1.9])   

% Plot Expected Arbitrage Benefit
figure
hold on
plot(d,f)
plot(d,0*f,'k--','LineWidth',2)
title('Expected Arbitrage Benefit')
xlabel('Dividend')
ylabel('Profit')
xlim([0.1 1.9])   

% Plot Response Residual
figure
hold on
plot(d,resid)
plot(d,0*resid,'k--','LineWidth',2)
title('Response Residual')
xlabel('Dividend')
ylabel('Residual')
xlim([0.1 1.9])   


%% SIMULATION

% Intizlizations
nper = 50;
nrep = 10000;
dsim = zeros(nrep,nper);
psim = zeros(nrep,nper);
dinit = ones(nrep,1);

% Simulate Model Not Using remsimul
rand('seed',0.945);
dd = dinit;
for ip=1:nper+1
  pp = interp1(d,p,dd,'cubic');
  dsim(:,ip,:) = dd;
  psim(:,ip,:) = pp;
  if ip<nper+1
    dd = dbar+gamma*(dd-dbar)+e(discrand(nrep,w),:);
  end
end

% Simulate Model Using remsimul
rand('seed',0.945);
[dsim,psim] = remsimul(model,basis,nper,dinit,d,p);

% Plot Simulated and Expected State
figure
plot(0:nper,dsim(1:3,:),0:nper,mean(dsim),'k')
title('Simulated and Expected Asset Dividend')
xlabel('Period')
ylabel('Dividend')

% Plot Simulated and Expected Response
figure
plot(0:nper,psim(1:3,:),0:nper,mean(psim),'k')
title('Simulated and Expected Asset Price')
xlabel('Period')
ylabel('Price')

% Compute Ergodic Moments
davg = mean(dsim(:)); 
pavg = mean(psim(:)); 
dstd = std(dsim(:)); 
pstd = std(psim(:)); 

% Print Ergodic Means and Standard Deviations
fprintf('Ergodic Means\n') 
fprintf('   Asset Dividend = %5.2f\n'    ,davg)
fprintf('   Asset Price    = %5.2f\n\n'  ,pavg)
fprintf('Ergodic Standard Deviations\n') 
fprintf('   Asset Dividend = %5.2f\n'    ,dstd)
fprintf('   Asset Price    = %5.2f\n\n'  ,pstd)

% Compute and Plot Ergodic Distribution
[qq,dd] = ksdensity(dsim(:),'support','positive','bandwidth',0.04);
figure
plot(dd,qq)
title('Ergodic Dividend Distribution')
xlabel('Dividend')
ylabel('Probability')
xlim([0 2])   

% Compute and Plot Ergodic Response Distribution
[qq,pp] = ksdensity(psim(:),'support','positive','bandwidth',0.04);
figure
plot(pp,qq)
title('Ergodic Price Distribution')
xlabel('Price')
ylabel('Probability')
xlim([5 15])  


%% Save Plots as EPS Files
printfigures(mfilename,7)


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

function [out1,out2,out3,out4] = func(flag,d,p,dn,pn,e,delta,beta,dbar,gamma)

ns = length(d);
switch flag
  case 'b'
    out1 = zeros(ns,1)-inf;
    out2 = zeros(ns,1)+inf;
  case 'f'
    u = d.^(-beta);
    un = dn.^(-beta);
    out1 = p.*u-delta*(pn+dn).*un;
    out2 = u;
    out3 = beta*delta*(pn+dn).*un./dn-delta*un;
    out4 = -delta*un;
  case 'g'
    out1 = dbar+gamma*(d-dbar)+e;
    out2 = zeros(ns,1);
end
%%