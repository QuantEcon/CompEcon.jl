function demdp19alt

%% DEMDP19ALT Credit With Strategic Default

% An infinitely-lived agent subject to employment uncertatintly must decide
% how much to consume, borrow, and save and whether to default on debt
% obligations, At the beginning of any period, the agent is either employed
% or unemployed and either ``credit-worthy'' or ``credit-unworthy''. The
% agent may save an unlimited amount and, if credit-worthy, may borrow up
% to a specified amount. A credit-worthy agent who is carrying debt may
% chose to default on her debt obligations. A defaulting credit-worthy
% agent completely erases her debt, but is immediately declared
% credit-unworthy and is barred from borrowing again until her credit is
% reinstated, which may occur from one period to the next a specified
% probability.
%
% This demo solves the collocation equation via Newton's method using
% dpsolve. The accompanying demo demdp19 solves the collcation equation
% directly using Broyden's method.
%
% States
%      s        assets (s<0 indicates debt)
%      i1       employment state (1=unemployed, 2=employed)
%      i2       credit state (1=unworthy, 2=worthy)
% Compound discrete state index by
%      i=1      same as i1=1 & i2=1 Unemployed unworthy agent
%      i=2      same as i1=2 & i2=1 Employed unworthy agent
%      i=3      same as i1=1 & i2=2 Unemployed worthy agent
%      i=4      same as i1=2 & i2=2 Employed worthy agent
% Actions
%      x        saving (x<0 indicates borrowing}
%      j        default decision (1=default, 2=not default)

% Preliminary tasks
% demosetup(mfilename)


%% FORMULATION
    
% Model Parameters
delta = 0.90;        % discount factor
alpha = 2.0;         % relative risk aversion
rb    = 0.10;        % interest rate on borrowing
rs    = 0.05;        % interest rate on saving
blim  = 0.5;         % borrowing limit
slim  = 1.0;         % saving limit
gamma = 0.07;        % unemployment rate
eta   = 2.0;         % expected duration of unemployment
y     = [0.5;1.0];   % income per employment state
sigma = 0.08;        % stigma penalty
mu    = 0.05;        % credit reinstatement probability

% Check Borrowing Limit
if min(y)<=rb*blim, warning('Inadmissible borrowing limit.'), end

% State Indices
[i1,i2] = indices(2,2);

% Utility Function
util = @(c) c.^(1-alpha)/(1-alpha);

% Employment State Transition Probabilities 
q      = zeros(2,2);
q(1,1) = 1-1/eta;
q(2,1) = gamma*(1-q(1,1))/(1-gamma);
q(1,2) = 1-q(1,1);
q(2,2) = 1-q(2,1);

% Compute Expected Lifetime Income
e = (eye(2,2)-delta*q)\y;

% Extend Income to Compound Discrete States
y = y(i1);

% Employment State Stationary Disbritution and Expected Visit Durations
p   = markov(q);
eta = 1./(1-diag(q));

% Print Employment State Transition Probabilities, Stationary Disbritution and Expected Visit Durations
fprintf('\n\n')
fprintf('Employment State Transition and Stationary Probabilities and Expected State Durations\n')
fprintf('           Unemployed   Employed   Stationary\n')
fprintf('Unemployed    %7.3f  %9.3f  %11.3f\n',q(1,:)',p(1))
fprintf('Employed      %7.3f  %9.3f  %11.3f\n',q(2,:)',p(2))
fprintf('Duration      %7.3f  %9.3f\n',eta(:))

% Action Contingent Probability Transition Matrix
qtilde = zeros(2,2,2,2,2);
qtilde(:,1,:,1,1) = q;
qtilde(:,1,:,2,1) = q*0;
qtilde(:,2,:,1,1) = q;
qtilde(:,2,:,2,1) = q*0;
qtilde(:,1,:,1,2) = q*(1-mu);
qtilde(:,1,:,2,2) = q*mu;
qtilde(:,2,:,1,2) = q*0;
qtilde(:,2,:,2,2) = q;
qtilde = reshape(qtilde,4,4,2);

% Discretize Continuous Action
m = 601;                                % number of continuous action nodes
X = nodeunif(m,-blim,slim);             % continuous action nodes
X = unique([X;0]);                      % add 0 anction node
m = length(X);                          % adjust number of nodes

% Approximation Structure
n     = 100;                            % number of collocation nodes
smin  = -(1+rb)*blim;                   % minimum net assets
smax  =  (1+rs)*slim;                   % maximum net assets
basis = fundefn('lin',n,smin,smax);     % basis functions

% Model Structure
model.func = @func;                     % model functions
model.discount = delta;                 % discount factor
model.ds = 1;                           % dimension of continuous state
model.dx = 1;                           % dimension of continuous action
model.ni = 4;                           % number of discrete states
model.nj = 2;                           % number of discrete actions
model.q  = qtilde;                      % discrete state transition probabilities
model.X  = X;                           % discretize continuous action
model.params = {util,rb,rs,blim,slim,y,i2,sigma};	% function parameters


%% SOLUTION

% Solve Bellman Equation
[c,s,v,x,resid] = dpsolve(model,basis);

% Express Redidual as Proportion of Expected Lifetime Income
resid = resid./(funeval(c,basis,s,1)*diag([e;e]));

% Compute Unconditional Saving Policy
xopt = x(:,:,1);
for i=1:4
  if i2(i)==2
    ind = find(v(:,i,2)>v(:,i,1));
    xopt(ind,i) = x(ind,i,2);
  end
end
n = length(s);

% Compute Asset Level at Which Agent Defaults
sdef = zeros(4,1)-inf;
for i=3:4
  if v(1,i,1)>v(1,i,2)
    ind = find(s<=0);
    sdef(i) = interp1(v(ind,i,2)-v(ind,i,1),s(ind),0,'linear','extrap');
  end
end

% Compute Asset Level at Which Agent Maximizes Debt
smaxdebt = [max([-inf;s(xopt(:,1)<=0)]); ...
            max([-inf;s(xopt(:,2)<=0)]); ...
            max([-inf;s(xopt(:,3)<=-blim)]); ...
            max([-inf;s(xopt(:,4)<=-blim)])];
smaxdebt = reshape(smaxdebt,2,2);

% Print Default and Max-out Levels
fprintf('Default and Max-out Asset Levels\n')
fprintf('                     Unemployed   Employed\n')
fprintf('Default                %7.3f  %9.3f\n',sdef(3:4))
fprintf('Max-out - Unworthy     %7.3f  %9.3f\n',smaxdebt(:,1))
fprintf('Max-out - Worthy       %7.3f  %9.3f\n',smaxdebt(:,2))
fprintf('\n')

% Plot Optimal Policy - Unworthy Agent
figure
hold on
plot(s,xopt(:,1:2),'LineWidth',4)
plot(smaxdebt(:,1),0,'*','LineWidth',6)
xlim([0 smax])
legend('Unemployed','Employed','Location','N')
legend('boxoff')
title('Optimal Saving - Unworthy Agent')
xlabel('Assets')
ylabel('Saving')

% Plot Optimal Policy - Worthy Agent
figure
hold on
plot(s,xopt(:,3:4),'LineWidth',4)
plot(smaxdebt(:,2),-blim,'*','LineWidth',6)
plot(sdef(3:4),-blim,'*','LineWidth',6)
xlim([smin smax])
legend('Unemployed','Employed','Location','N')
legend('boxoff')
title('Optimal Saving - Worthy Agent')
xlabel('Assets')
ylabel('Saving')

% Plot Considitional Value Functions for Unemployed Worthy Agent
figure
plot(s,squeeze(v(:,3,:)),'LineWidth',4)
xlim([smin smax])
legend('Default','Do Not Default','Location','Best')
legend('boxoff')
title('Conditional Value Functions - Unemployed Worthy Agent')
xlabel('Assets')
ylabel('Value')

% Plot Considitional Value Functions for Employed Worthy Agent
figure
plot(s,squeeze(v(:,4,:)),'LineWidth',4)
xlim([smin smax])
legend('Default','Do Not Default','Location','Best')
legend('boxoff')
title('Conditional Value Functions - Employed Worthy Agent')
xlabel('Assets')
ylabel('Value')

% Plot Value Function
figure
% Compute Unconditional Value Function
vopt = max(v(:,:,1),v(:,:,2));
plot(s,vopt,'LineWidth',4)
xlim([smin smax])
legend('Unworthy-Unemployed','Unworthy-Employed','Worthy-Unemployed','Worthy-Employed','Location','SE')
legend('boxoff')
title('Value Function')
xlabel('Assets')
ylabel('Value')

% Plot Residual
figure
plot(s,resid,'LineWidth',4)
xlim([smin smax])
legend('Unworthy-Unemployed','Unworthy-Employed','Worthy-Unemployed','Worthy-Employed','Location','SE')
legend('boxoff')
title('Bellman Equation Residual')
xlabel('Assets')
ylabel('Proportion of Expected Lifetime Income')


%% SIMULATION

% Compute Asset Transitions
snext = (1+rb)*min(xopt,0) + (1+rs)*max(xopt,0);
snext = reshape(snext,n,2,2);

%  Simulation Control Parameters
nrep = 4000;           	% number of agents simulated
nper =  500;           	% number of periods simulated
nwrm =   50;          	% number of `warm-up' periods
rng('default')        	% seed random number generator

% Initialize History Arrays
i1simhist = zeros(nrep,nper);
i2simhist = zeros(nrep,nper);
ssimhist  = zeros(nrep,nper);
xsimhist  = zeros(nrep,nper);
jsimhist  = zeros(nrep,nper);

% Intialize: all agents employed, credit-worthy, no assets
i1sim = 2*ones(nrep,1);
i2sim = 2*ones(nrep,1);
ssim  = zeros(nrep,1);

% Simulate
for t=1:nper
  % determine default decision
  jsim  = 2*ones(nrep,1);
  for i1=1:2
    jsim(i1sim==i1&i2sim==2&ssim<sdef(2+i1)) = 1;
  end
  % store states and action in history arrays
  i1simhist(:,t) = i1sim;
  i2simhist(:,t) = i2sim;
  ssimhist(:,t)  = ssim;
  jsimhist(:,t)  = jsim;
  if t<nper
    % update employment state
    i1simnew = markovsim(i1sim,q);
    % update creditworthiness state
    i2simnew = zeros(nrep,1);
    i2simnew(i2sim==1) = 1;
    i2simnew(i2sim==1&rand(nrep,1)<mu) = 2;
    i2simnew(i2sim==2) = jsim(i2sim==2);
    % update net assets
    ssimnew  = zeros(nrep,1);
    for i1=1:2
      for i2=1:2
        ind = find(i1sim==i1&i2sim==i2);
        if ~isempty(ind)
          ssimnew(ind) = interp1(s,snext(:,i1,i2),ssim(ind));
        end
      end
    end
    % change in assets
    xsimhist(:,t) = ssimnew-ssim;
    % update states
    i1sim = i1simnew;
    i2sim = i2simnew;
    ssim  = ssimnew;
  end
end

% Plot Simulated State Path
figure
hold on
plot(0:40,mean(ssimhist(:,1:41)),'k')
title('Expected Assets')
xlabel('Period')
ylabel('Assets')

% Plot Simulated State Path
figure
hold on
plot(0:40,2-mean(i2simhist(:,1:41)),'k')
title('Expected Delinquency Rate')
xlabel('Period')
ylabel('Percent')

% Post-Warmup Period Observations to Approximate Ergodic Distribution
i1simhist = i1simhist(:,nwrm:nper); i1simhist = i1simhist(:); 
i2simhist = i2simhist(:,nwrm:nper); i2simhist = i2simhist(:); 
ssimhist  = ssimhist(:,nwrm:nper);  ssimhist  = ssimhist(:);  
xsimhist  = xsimhist(:,nwrm:nper);  xsimhist  = xsimhist(:);
jsimhist  = jsimhist(:,nwrm:nper);  jsimhist  = jsimhist(:);  

% Compute Ergodic Means and Discrete State Visit Probabilities
savg = zeros(4,1);
xavg = zeros(4,1);
prop = zeros(4,1);
savg(1) = mean(ssimhist(i1simhist==1));
xavg(1) = mean(xsimhist(i1simhist==1));
prop(1) = mean(i1simhist==1);
savg(2) = mean(ssimhist(i1simhist==2));
xavg(2) = mean(xsimhist(i1simhist==2));
prop(2) = mean(i1simhist==2);
savg(3) = mean(ssimhist(i2simhist==1));
xavg(3) = mean(xsimhist(i2simhist==1));
prop(3) = mean(i2simhist==1);
savg(4) = mean(ssimhist(i2simhist==2));
xavg(4) = mean(xsimhist(i2simhist==2));
prop(4) = mean(i2simhist==2);
savgall = mean(ssimhist);
xavgall = mean(xsimhist);

% Compute Ergodic Standard Deviations
sstd = zeros(4,1);
xstd = zeros(4,1);
sstd(1) = std(ssimhist(i1simhist==1));
xstd(1) = std(xsimhist(i1simhist==1));
sstd(2) = std(ssimhist(i1simhist==2));
xstd(2) = std(xsimhist(i1simhist==2));
sstd(3) = std(ssimhist(i2simhist==1));
xstd(3) = std(xsimhist(i2simhist==1));
sstd(4) = std(ssimhist(i2simhist==2));
xstd(4) = std(xsimhist(i2simhist==2));
sstdall = std(ssimhist);
xstdall = std(xsimhist);

% Compute Ergodic Default Rates
rdef = zeros(4,1);
rdef(1) = sum(i1simhist==1&ssimhist<0&jsimhist==1)/sum(i1simhist==1&ssimhist<0);
rdef(2) = sum(i1simhist==2&ssimhist<0&jsimhist==1)/sum(i1simhist==2&ssimhist<0);
rdef(3) = sum(i2simhist==1&ssimhist<0&jsimhist==1)/sum(i2simhist==1&ssimhist<0);
rdef(4) = sum(i2simhist==2&ssimhist<0&jsimhist==1)/sum(i2simhist==2&ssimhist<0);
rdefall = sum(i2simhist==2&ssimhist<0&jsimhist==1)/sum(i2simhist==2&ssimhist<0);

% Print Ergodic Means, Default Rates and Visit Probabilities
fprintf('\nErgodic Means, Default Rates and Visit Probabilities\n')
fprintf('           Unemployed   Employed   Unworthy     Worthy       All\n')
fprintf('Assets        %7.3f  %9.3f  %9.3f  %9.3f  %9.3f\n',savg(:),savgall)
fprintf('Saving        %7.3f  %9.3f  %9.3f  %9.3f  %9.3f\n',xavg(:),xavgall)
fprintf('Default Rate  %7.3f  %9.3f  %9.3f  %9.3f  %9.3f\n',rdef(:),rdefall)
fprintf('Probability   %7.3f  %9.3f  %9.3f  %9.3f  %9.3f\n',prop(:),sum(prop(:))/2)

% Print Ergodic Standard Deviations
fprintf('\nErgodic Standard Deviations\n')
fprintf('           Unemployed   Employed   Unworthy     Worthy       All\n')
fprintf('Assets        %7.3f  %9.3f  %9.3f  %9.3f  %9.3f\n',sstd(:),sstdall)
fprintf('Saving        %7.3f  %9.3f  %9.3f  %9.3f  %9.3f\n',xstd(:),xstdall)

% % Compute and Plot Ergodic Wealth Distribution
% [qq,ss] = ksdensity(ssimhist+y(i1simhist),'bandwidth',0.01);
% figure
% plot(ss,qq)
% title('Ergodic Wealth Distribution')
% xlabel('Wealth')
% ylabel('Probability')


%% Save Plots as EPS Files
% printfigures(mfilename,8)


%% DPSOLVE FUNCTION FILE
%
%    User-supplied function called by dpsolve that returns the bound,
%    reward, and continuous state transition function values and
%    derivatives with respect to the continuous action x at an arbitrary
%    number ns of states and actions aaccording to the format
%      [out1,out2,out3] = func(flag,s,x,i,j,in,e,<params>)
%    For further information regarding input format, type > help dpsolve or
%    reference the DPSOLVE Handbook.

function [out1,out2] = func(flag,s,x,i,j,in,e,util,rb,rs,blim,slim,y,i2,sigma)
i2 = i2(i)-1;
j  = j-1;
if j==0
  s = max(s,0);
end
n = length(s);
switch flag
  case 'b'
    if j==0||i2==0
      out1 = zeros(n,1);
    else
      out1 = zeros(n,1)-blim;
    end
    out2 = zeros(n,1)+slim;
  case 'f'
    out1 = -inf+zeros(n,1);
    c = s + y(i) - x;
    if j==0||i2==0
      k = find(c>0&x>=0);
    else
      k = find(c>0);
    end
    out1(k) = util(c(k))-(1-i2)*sigma;
  case 'g'
    out1 = (1+rb)*min(x,0) + (1+rs)*max(x,0);
end
%%