function demdp15

%% DEMDP15 Saving with Transactions Costs

% An infinitely-lived agent subject to income uncertainty must decide how
% much to consume, borrow, and save, given that he incurs a fixed
% transaction cost whenever he adjusts his assets and given interest rates
% on debt and deposits differ.

% States
%      s        beginning assets (s<0 indicates debt)
%      i        income state
% Actions
%      send     ending assets
%      j        transaction decision (no=1, yes=2)

% Preliminary tasks
demosetup(mfilename)

   
%% FORMULATION

% Model Parameters
delta = 0.95;        % discount factor
alpha = 2.0;         % relative risk aversion
rb    = 0.15;        % interest rate on debt
rs    = 0.05;        % interest rate on deposits
tau   = 0.03;        % fixed transactions cost
smin  = -0.5;        % minimum assets (debt limit)
smax  =  1.0;        % maximum assets (deposit limit)
y     = [1.0;0.6];   % high-low incomes
w     = [0.7;0.3];   % probabilities of high-low incomes
ny    = length(y);   % number of income states

% Check Borrowing Limit
if min(y)+rb*smin<=0, warning('Inadmissible borrowing limit.'), end

% Utility Function
util = @(c) c.^(1-alpha)/(1-alpha);

% Discretize Admissible Asset Holdings
n = 601;
s = nodeunif(n,smin,smax);
s = unique([s;0]);
n = length(s);

% State Matrix for Discrete Optimization
S = s(:,ones(1,n));

% Net Interest
Is = rb*min(s,0) + rs*max(s,0);
IS = rb*min(S,0) + rs*max(S,0);


%% SOLUTION

% Intialize
v  = zeros(n,ny);
v1 = zeros(n,ny);
v2 = zeros(n,ny);

% Perform some Function Iterations on Collocation Equations
for it=1:100
  vold = v;
  Evnext1 = v*w;
  Evnext2 = Evnext1(:,ones(1,n))';
  for i=1:ny
    % do not transact
    v1(:,i) = util(y(i)+Is) + delta*Evnext1;
    % transact
    C = y(i) + IS + S - S' - tau;
    U = util(C);
    U(C<=0) = -inf;
    v2(:,i) = max(U+delta*Evnext2,[],2);
  end
  v = max(v1,v2);
  change = norm(v(:)-vold(:));
  fprintf ('%4i %10.1e\n',it,change)
  if change<1.e-8, break, end
end

% Solve Collocation Equation by Broyden's Method
tic
optset('broyden','initi',1)
optset('broyden','showiters',1)
optset('broyden','tol',1e-10)
optset('broyden','maxit',200)
optset('broyden','maxsteps',5)
v = broyden(@resid,v(:),n,ny,util,s,S,Is,IS,y,w,tau,delta);
toc

% Recover Optimal Policy
[r,send] = resid(v,n,ny,util,s,S,Is,IS,y,w,tau,delta);

% Plot Optimal Policy
figure
hold on
title('Optimal Saving Policy')
plot(s,send,'LineWidth',4)
legend('High Income','Low Income','Location','Best') 
legend('boxoff')
plot(s,s,'k:','LineWidth',1)
plot(s,0,'k:','LineWidth',2)
xlim([smin smax])
xlabel('Beginning Assets')
ylabel('Ending Assets')

% Plot Value Function
figure
hold on
title('Value Function')
v = reshape(v,n,ny);
plot(s,v,'LineWidth',4)
legend('High Income','Low Income','Location','Best') 
legend('boxoff')
xlim([smin smax])
xlabel('Beginning Assets')
ylabel('Value')


%% SIMULATION

%  Simulation Control Parameters
nrep = 4000;           	% number of agents simulated
nper = 500;           	% number of periods simulated
nwrm = 100;          	% number of `warm-up' periods
rng('default')        	% seed random number generator

% Initialize History Arrays
isimhist = zeros(nrep,nper);   % income state
ssimhist = zeros(nrep,nper);   % beginning assets
xsimhist = zeros(nrep,nper);   % change in assets

% Intialize: all agents employed, no assets
isim = 2*ones(nrep,1);
ssim = zeros(nrep,1);

% Simulate
for t=1:nper
  % store states in history arrays
  isimhist(:,t) = isim;
  ssimhist(:,t) = ssim;
  if t<nper
    % new employment state
    isimnew = discrand(nrep,w);
    % new ending assets
    ssimnew  = zeros(nrep,1);
    for i=1:ny
      ind = find(isim==i);
      if ~isempty(ind)
        ssimnew(ind) = interp1(s,send(:,i),ssim(ind));
      end
    end
    % change in assets
    xsimhist(:,t) = ssimnew-ssim;
    % update states
    isim = isimnew;
    ssim = ssimnew;
  end
end

% Plot Simulated State Path
figure
hold on
T = 30;
plot(0:T,ssimhist(1:3,1:T+1),0:T,mean(ssimhist(:,1:T+1)),'k')
title('Simulated and Expected Beginning Assets')
xlabel('Period')
ylabel('Assets')

% Retain Post-Warmup Period Observations 
isimhist = isimhist(:,nwrm:nper); isimhist = isimhist(:);
ssimhist = ssimhist(:,nwrm:nper); ssimhist = ssimhist(:);
xsimhist = xsimhist(:,nwrm:nper); xsimhist = xsimhist(:);

% Compute Ergodic Mean Beginning Assets, Change in Assets and Visit Probabilities, per Income State
savg = zeros(ny,1);
xavg = zeros(ny,1);
prop = zeros(ny,1);
for i=1:ny
  savg(i) = mean(ssimhist(isimhist==i));
  xavg(i) = mean(xsimhist(isimhist==i));
  prop(i) = mean(isimhist==i);
end
xavgall = mean(xsimhist);
savgall = mean(ssimhist);

% Compute Ergodic Standard Deviations of Assets and Asset Changes, per Income State
sstd = zeros(ny,1);
xstd = zeros(ny,1);
for i=1:ny
  sstd(i) = std(ssimhist(isimhist==i)+xsimhist(isimhist==i));
  xstd(i) = std(xsimhist(isimhist==i));
end
sstdall = std(ssimhist+xsimhist);
xstdall = std(xsimhist);

% Print Ergodic Mean Assets and Visit Probabilities, per Income State
fprintf('\nErgodic Mean Assets and Visit Probabilities, per Income State\n')
fprintf('                     High       Low\n')
fprintf('                   Income     Income       All\n')
fprintf('Ending Assets     %7.3f  %9.3f %9.3f\n',savg+xavg,savgall+xavgall)
fprintf('Beginning Assets  %7.3f  %9.3f %9.3f\n',savg,savgall)
fprintf('Change in Assets  %7.3f  %9.3f %9.3f\n',xavg,xavgall)
fprintf('Probability       %7.3f  %9.3f %9.3f\n',prop,sum(prop))

% Print Ergodic Standard Deviation of Assets, per Income State
fprintf('\nErgodic Standard Devitations of Assets, per Income State\n')
fprintf('                     High       Low\n')
fprintf('                   Income     Income       All\n')
fprintf('Beginning Assets  %7.3f  %9.3f %9.3f\n',sstd,sstdall)
fprintf('Change in Assets  %7.3f  %9.3f %9.3f\n',xstd,xstdall)


%% Save Plots as EPS Files
if tau==0&&~(rb==rs)
  printfigures('demdp15not',3)
elseif ~(tau==0)&&rb==rs
  printfigures('demdp15eqr',3)
elseif tau==0&&rb==rs
  printfigures('demdp15noteqr',3)
else
  printfigures(mfilename,3)
end


%% RESIDUAL FUNCTION

function [r,send] = resid(v,n,ny,util,s,S,Is,IS,y,w,tau,delta)

% Reshape and Intialize
v  = reshape(v,n,ny);
v1 = zeros(n,ny);
v2 = zeros(n,ny);
send = zeros(n,ny);

% Evaluate RHS of Bellman Equation
Evnext1 = v*w;
Evnext2 = Evnext1(:,ones(1,n))';
for i=1:ny
  % do not transact
  v1(:,i) = util(y(i)+Is) + delta*Evnext1;
  % transact
  C = y(i) + IS + S - S' - tau;
  U = util(C);
  U(C<=0) = -inf;
  [v2(:,i),isn] = max(U+delta*Evnext2,[],2);
  send(:,i) = s(isn);
end
r = v(:) - max(v1(:),v2(:));

for i=1:ny
  ind = find(v1(:,i)>v2(:,i));
  send(ind,i) = s(ind);
end
%%