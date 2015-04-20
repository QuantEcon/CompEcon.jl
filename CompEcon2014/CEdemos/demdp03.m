function demdp03

%% DEMDP03 Industry Entry-Exit Model
%
% A profit maximizing firm must decide whether to operate or shut down, 
% given its short-run profitability, subject to transactions costs.
%
% States
%     p       current short-run profit
%     i       active (1) or idle (2) last period
% Actions
%     j       active (1) or idle (2) this period
% Parameters
%     pbar    long-run mean profit
%     gamma   profit autoregressive coefficient
%     kappa   cost of reopenning idle firm
%     sigma   standard deviation of profit shock
%     delta   discount factor
  
% Preliminary tasks
demosetup(mfilename)


%% FORMULATION
  
% Model Parameters
pbar   = 1.0;                               % long-run mean profit
gamma  = 0.7;                               % profit autoregressive coefficient
kappa  =  10;                             	% cost of reopenning idle firm
sigma  = 1.0;                               % standard deviation of profit shock
delta  = 0.9;                             	% discount factor
  
% Continuous State Shock Distribution
m = 5;                                      % number of profit shocks
[e,w] = qnwnorm(m,0,sigma.^2);              % profit shocks and probabilities

% Approximation Structure
n = 250;                                    % number of collocation nodes
pmin = -20;                                 % minimum profit
pmax =  20;                                 % maximum profit
basis = fundefn('spli',n,pmin,pmax);        % basis functions
  
% Model Structure
model.func = @func;                         % model functions
model.params = {pbar gamma kappa};          % function parameters
model.discount = delta;                     % discount factor
model.ds = 1;                               % dimension of continuous state
model.dx = 0;                               % dimension of continuous action
model.ni = 2;                               % number of discrete states
model.nj = 2;                               % number of discrete actions
model.e  = e;                              	% unit profit contribution shocks
model.w  = w;                              	% unit profit shock probabilities
model.h  = [1 1;2 2];                       % deterministic discrete state transitions


%% SOLUTION

% Solve Bellman Equation
[c,sr,vr,xr,resid] = dpsolve(model,basis); 

% Plot Action-Contingent Value Functions
figure
hold on
plot(sr,[vr(:,1,1) vr(:,2,1) vr(:,2,2)])
legend('Keep Active Firm Open','Reopen Idle Firm','Shut Down','location','NW')
title('Action-Contingent Value Functions')
xlabel('Potential Profit')
ylabel('Value of Firm')

%  Compute and Plot Critical Profit
for i=1:2
  pcrit = interp1(vr(:,i,1)-vr(:,i,2),sr,0);
  vcrit  = interp1(sr,vr(:,i,1),pcrit);
  plot([pcrit pcrit],[-40 vcrit],'k--')
  bullet(pcrit,-40)
  bullet(pcrit,vcrit)
  ftext(pcrit,-40,['\bf p^*_' int2str(2-i)],'left','bottom')
  if i==1
    fprintf('Profit Exit  = %5.2f\n'    ,pcrit)  
  else
    fprintf('Profit Entry = %5.2f\n\n'  ,pcrit)
  end  
end

% Plot Residual
figure
hold on
plot(sr,100*resid./max(vr,[],3))
plot(sr,0*sr,'k--')
% ylim([-0.16 0.16])
legend('Active','Idle','location','SE')
title('Bellman Equation Residual')
xlabel('Potential Profit')
ylabel('Percent Residual')


%% SIMULATION

% Simulate Model
rand('seed',0.945);
nper = 51; 
nrep = 50000;
pinit = pbar*ones(nrep,1);
iinit = ones(nrep,1);
[ssim,xsim,isim] = dpsimul(model,basis,nper,pinit,iinit,sr,vr);

% Compute Ergodic Moments
isim = 2 - isim; % Convert to 0=idle, 1=active
savg = mean(ssim(:,nper));
iavg = mean(isim(:,nper));
sstd = std(ssim(:,nper));
istd = std(isim(:,nper));

% Print Ergodic Moments
fprintf('Ergodic Means\n') 
fprintf('   Profit Contribution  = %5.2f\n'    ,savg)
fprintf('   Activity             = %5.2f\n\n'  ,iavg)
fprintf('Ergodic Standard Deviations\n') 
fprintf('   Profit Contribution  = %5.2f\n'    ,sstd)
fprintf('   Activity             = %5.2f\n\n'  ,istd)

% Plot Simulated and Expected Continuous State Path
figure
hold on
plot(0:nper-1,ssim(1:3,:))
plot(0:nper-1,mean(ssim),'k')
title('Simulated and Expected Profit Contribution')
xlabel('Period')
ylabel('Profit Contribution')

% Plot Expected Discrete State Path
figure
hold on
% plot(0:nper,isim(1:3,:))
plot(0:nper-1,mean(isim),'k')
title('Probability of Operation')
xlabel('Period')
ylabel('Probability')


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

function out = func(flag,p,x,i,j,in,e,pbar,gamma,kappa)

switch flag
case 'f'
  out = p.*(j==1) - kappa.*(i==2).*(j==1);
case 'g'
  out = pbar+gamma*(p-pbar)+e;
end
%%