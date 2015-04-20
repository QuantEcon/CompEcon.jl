function demdp02

%% DEMDP02 Asset Replacement Model
%
% Profit-maximizing entrepreneur must decide when to replace an aging asset.
%
% States
%     p       unit profit
%     a       asset age (1..A)
% Actions
%     j       keep(1) or replace(2) asset
% Parameters
%     A       maximum asset age 
%     alpha   production function coefficients
%     kappa   net replacement cost
%     pbar    long-run mean unit profit
%     gamma   unit profit autoregression coefficient
%     sigma   standard deviation of unit profit shock
%     delta   discount factor 

% Preliminary tasks
demosetup(mfilename)


%% FORMULATION
  
% Model Parameters
A       = 6;                                % maximum asset age 
alpha   = [50 -2.5 -2.5];                   % production function coefficients
kappa   = 40;                               % net replacement cost
pbar    = 1;                                % long-run mean unit profit
gamma   = 0.5;                              % unit profit autoregression coefficient
sigma   = 0.15;                             % standard deviation of unit profit shock
delta   = 0.9;                              % discount factor 

% Continuous State Shock Distribution
m = 5;                                      % number of unit profit shocks
[e,w] = qnwnorm(m,0,sigma^2);               % unit profit shocks and probabilities

% Deterministic Discrete State Transitions
h = [2:A 1; ones(1,A)];

% Approximation Structure
n  = 200;                                   % number of collocation nodes
pmin = 0;                                   % minimum unit profit
pmax = 2;                                   % maximum unit profit
basis = fundefn('spli',n,pmin,pmax);        % basis functions

% Model Structure
model.func = @func;                         % model functions
model.params = {A alpha kappa pbar gamma};  % function parameters
model.discount = delta;                   	% discount factor
model.ds = 1;                               % dimension of continuous state
model.dx = 0;                               % dimension of continuous action
model.ni = A;                               % number of discrete states
model.nj = 2;                               % number of discrete actions
model.e  = e;                            	% unit profit shocks
model.w  = w;                              	% unit profit shock probabilities
model.h  = h;                              	% deterministic discrete state transitions


%% SOLUTION

% Solve Bellman Equation
[c,pr,vr,xr,resid] = dpsolve(model,basis); 

% Plot Action-Contingent Value Functions
figure
hold on
plot(pr,squeeze(vr(:,:,1)))
legend('Keep a=1','Keep a=2','Keep a=3','Keep a=4','Keep a=5','Replace','Location','NW')
title('Action-Contingent Value Functions')
xlabel('Net Unit Profit')
ylabel('Value')

% Compute and Plot Critical Unit Profit Contributions
fprintf('Critical Replacement Profit\n')
for a=1:A-1
  pcrit = interp1(vr(:,a,1)-vr(:,a,2),pr,0);
  vcrit = interp1(pr,vr(:,a,1),pcrit);
  if isnan(pcrit), continue, end
  bullet(pcrit,150)
  bullet(pcrit,vcrit)
  plot([pcrit pcrit],[150 vcrit],'k--')
  ftext(pcrit,150,['\bf p^*_' int2str(a)],'left','bottom')
  fprintf('   Age %2i  Profit %5.2f\n'  ,a,pcrit)
end
fprintf('\n')

% Plot Residual
figure
hold on
plot(pr,100*resid./max(vr,[],3))
plot(pr,0*pr,'k--')
legend('a=1','a=2','a=3','a=4','a=5','a=6','Location','NE')
title('Bellman Equation Residual')
xlabel('Net Unit Profit')
ylabel('Percent Residual')


%% SIMULATION

% Simulate Model
rand('seed',0.945);
nper = 51; 
nrep = 10000;
sinit = pbar*ones(nrep,1);
iinit = ones(nrep,1);
[ssim,xsim,isim,jsim] = dpsimul(model,basis,nper,sinit,iinit,pr,vr);

% Compute Ergodic Moments
savg = mean(ssim(:));
iavg = mean(isim(:));
sstd = std(ssim(:));
istd = std(isim(:));

% Print Ergodic Moments
fprintf('Ergodic Means\n') 
fprintf('   Price        = %5.2f\n'      ,savg)
fprintf('   Age          = %5.2f\n\n'    ,iavg)
fprintf('Ergodic Standard Deviations\n') 
fprintf('   Price        = %5.2f\n'      ,sstd)
fprintf('   Age          = %5.2f\n\n'    ,istd)

% Plot Simulated and Expected Continuous State Path
figure
hold on
plot(0:nper-1,ssim(1:3,:))
plot(0:nper-1,mean(ssim),'k')
title('Simulated and Expected Price')
xlabel('Period')
xlabel('Net Unit Profit')

% Plot Expected Discrete State Path
figure
hold on
plot(0:nper-1,mean(isim),'k')
title('Expected Machine Age')
xlabel('Period')
ylabel('Age')


%% Save Plots as EPS Files
printfigures(mfilename,4)


%% DPSOLVE FUNCTION FILE
%
%    User-supplied function called by dpsolve that returns the bound,
%    reward, and continuous state transition function values and
%    derivatives with respect to the continuous action x at an arbitrary
%    number ns of states and actions aaccording to the format
%      [out1,out2,out3] = func(flag,s,x,i,j,in,e,<params>)
%    For further information regarding input format, type > help dpsolve or
%    reference the DPSOLVE Handbook.

function out = func(flag,p,x,a,j,in,e,A,alpha,kappa,pbar,gamma)

switch flag
case 'f'
  if j==2||a==A
    out = p*50-kappa;
  else
    out = p*(alpha(1)+alpha(2)*a+alpha(3)*a.^2);
  end
case 'g'
  out = pbar + gamma*(p-pbar) + e;
end
%%