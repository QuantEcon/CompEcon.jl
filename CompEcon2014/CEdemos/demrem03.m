function demrem03

%% DEMREM03 Government Price Support Model
%
% Solves standard government commodity price support model.
%
% States
%     s         supply
% Responses
%     a         acreage planted
%     z         government stocks
% Parameters
%     zbar      maximum stocks
%     pbar      government support price
%     gamma     inverse consumption demand elasticity
%     beta      acreage supply elasticity
%     yvol      yield volatility
%     delta     discount factor

% Preliminary tasks
demosetup(mfilename)


%% FORMULATION

% Model Parameters
zbar  = 1.0;                                    % maximum stocks
pbar  = 0.9;                                    % government support price
gamma = 2.0;                                    % inverse consumption demand elasticity
beta  = 0.5;                                    % acreage supply elasticity
yvol  = 0.3;                                    % yield volatility
delta = 0.9;                                    % discount factor

% Shock Distribution
m = 11;                                         % number of yield shocks
[y,w] = qnwlogn(m,-0.5*yvol^2,yvol^2);          % log normal nodes and weights

% Approximation Structure
n      = 501;                                   % degree of approximation
smin   = 0.3;                                   % minimum supply
smax   = 3.5;                                   % maximum supply
basis  = fundefn('spli',n,smin,smax);           % basis functions
s = funnode(basis);

% Model Structure
model.func = @func;                             % model functions
model.params = {delta zbar pbar gamma beta};    % function parameters
model.ds = 1;                                   % dimension of state s
model.dx = 2;                                   % dimension of response x
model.e = y;                                    % shocks
model.w = w;                                    % probabilities


%% SOLUTION

% Consumption Demand at Support Price
dstar = pbar^(-1/gamma);

% Compute Equilibrium Government Stocks and Price
z = min(max(s-dstar,0),zbar);
p = (s-z).^(-gamma);

% Compute Equilibrium Acreage Planted Directly
maxit = 100;
a = ones(n,1);
for it=1:maxit
  Er = 0;
  Dr = 0;
  for j=1:m
    sn = z + a*y(j);
    zn = min(max(sn-dstar,0),zbar);
    pp = (sn-zn).^(-gamma);
    da = zeros(n,1);
    da(0>sn-dstar | sn-dstar>zbar) = y(j);
    dd = -gamma*((sn-zn).^(-gamma-1)).*da;
    Er = Er + w(j)*pp*y(j);
    Dr = Dr + w(j)*dd*y(j);
  end
  f = a - (delta*Er).^(1/beta);
  d = 1 - delta*(1/beta)*Dr.*(delta*Er).^((1/beta)-1);
  a = a - f./d;
  if norm(f)<1.e-10, break, end
end

% % Compute Rational Expectations Equilibrium Using remsolve
% x  = [a z];   
% [c,s,x,f] = remsolve(model,basis,x);
% a  = x(:,1);
% z  = x(:,2);
% p = (s-z).^(-gamma);

% Plot Equilibrium Response 1
figure
plot(s,a)
title('Equilibrium Acreage Planted')
xlabel('Supply')
ylabel('Acreage')
xlim([0 smax])

% Plot Equilibrium Response 2
figure
plot(s,z)
title('Equilibrium Government Stockholding')
xlabel('Supply')
ylabel('Stocks')
xlim([0 smax])  

% Plot Equilibrium Response 3
figure
hold on
plot(s,s.^(-gamma),'r')
plot(s,p)
legend('Total Demand','Consumption Demand')
legend boxoff
plot([0 smax],[pbar pbar],'k:')
title('Equilibrium Price')
xlabel('Quantity')
ylabel('Price')
xlim([0 smax])    
ylim([0 4])  


%% SIMULATION

% Initializations
nper  = 30;
nrep  = 10000;
sinit = ones(nrep,1);
ssim  = zeros(nrep,nper+1);
asim  = zeros(nrep,nper+1);
zsim  = zeros(nrep,nper+1);

% Simulate Model
rand('seed',0.945);
ss = sinit;
for ip=1:nper+1
  zz = min(max(ss-dstar,0),zbar);
  aa = interp1(s,a,ss,'cubic');
  ssim(:,ip,:) = ss;
  asim(:,ip,:) = aa;
  zsim(:,ip,:) = zz;
  if ip<nper+1
    ss = zz + aa.*y(discrand(nrep,w),:);
  end
end
[smin min(ss) max(ss) smax]
psim = (ssim-zsim).^(-gamma);

% Plot Simulated and Expected State 1
figure
plot(0:nper,ssim(1:2,:),0:nper,mean(ssim),'k')
title('Simulated and Expected Supply')
xlabel('Period')
ylabel('Supply')

% Plot Simulated and Expected Response 1
figure
plot(0:nper,asim(1:2,:),0:nper,mean(asim),'k')
title('Simulated and Expected Acreage Planted')
xlabel('Period')
ylabel('Acreage')

% Plot Simulated and Expected Response 2
figure
plot(0:nper,zsim(1:2,:),0:nper,mean(zsim),'k')
title('Simulated and Expected Government Stocks')
xlabel('Period')
ylabel('Stocks')

% Plot Simulated and Expected Response 3
figure
plot(0:nper,psim(1:2,:),0:nper,mean(psim),'k')
title('Simulated and Expected Price')
xlabel('Period')
ylabel('Price')

% Compute Ergodic Moments
savg = mean(ssim(:));
aavg = mean(asim(:));
zavg = mean(zsim(:));
pavg = mean(psim(:));
sstd = std(ssim(:));
astd = std(asim(:));
zstd = std(zsim(:));
pstd = std(psim(:));

% Print Ergodic Moments
fprintf('Ergodic Means\n') 
fprintf('   Supply   = %5.2f\n'    ,savg)
fprintf('   Acreage  = %5.2f\n'    ,aavg)
fprintf('   Stocks   = %5.2f\n'    ,zavg)
fprintf('   Price    = %5.2f\n\n'  ,pavg)
fprintf('Ergodic Standard Deviations\n') 
fprintf('   Supply   = %5.2f\n'    ,sstd)
fprintf('   Acreage  = %5.2f\n'    ,astd)
fprintf('   Stocks   = %5.2f\n'    ,zstd)
fprintf('   Price    = %5.2f\n\n'  ,pstd)

% Compute and Plot Ergodic State Distribution
[qq,ss] = ksdensity(ssim(:),'support','positive','bandwidth',0.14);
figure
plot(ss,qq)
title('Ergodic Supply Distribution')
xlabel('Supply')
ylabel('Probability')
xlim([0 4])  


%% Save Plots as EPS Files
printfigures(mfilename,8)


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

function [out1,out2,out3,out4] = func(flag,s,x,sn,xn,e,delta,zbar,pbar,gamma,beta)
a  = x(:,1);
z  = x(:,2);
ns = length(s);
ds = 1;
dx = 2;
switch flag
  case 'b'
    out1 = zeros(ns,2);                              
    out2 = [inf+zeros(ns,1) zbar+zeros(ns,1)];           
  case 'f'
    zn = xn(:,2);
    out1 = zeros(ns,dx);
    out2 = zeros(ns,dx,dx);
    out3 = zeros(ns,dx,ds);
    out4 = zeros(ns,dx,dx);
    p     = (s -z ).^(-gamma);
    pn    = (sn-zn).^(-gamma);
    pder  = -gamma*(s -z ).^(-gamma-1);
    pnder = -gamma*(sn-zn).^(-gamma-1);
    out1(:,1)   = delta*pn.*(sn-z)./a-a.^beta;
    out1(:,2)   = pbar-p;
    out2(:,1,1) = -delta*pn.*(sn-z)./(a.^2)-beta*a.^(beta-1);
    out2(:,1,2) = -delta*pn./a;
    out2(:,2,2) = pder;
    out3(:,1,1) = delta*(pn+pnder.*(sn-z))./a;
    out4(:,1,2) = -delta*pnder.*(sn-z)./a;
  case 'g'
    out2 = zeros(ns,ds,dx);
    out1 = z+a.*e;                            
    out2(:,1,1) = e;                                          
    out2(:,1,2) = ones(ns,1);                       
end
%%