function demdp04

%% DEMDP04 Job Search Model
%
% Infinitely-lived worker must decide whether to quit, if employed, or 
% search for a job, if unemployed, given prevailing market wages.
%
% States
%     w       prevailing wage
%     i       unemployed (1) or employed (2) at beginning of period
% Actions
%     j       idle (1) or active (i.e., work or search) (2) this period
% Parameters
%     v        benefit of pure leisure
%     wbar     long-run mean wage
%     gamma    wage reversion rate
%     p0       probability of finding job
%     p1       probability of keeping job
%     sigma    standard deviation of wage shock
%     delta    discount factor

% Preliminary tasks
demosetup(mfilename)


%% FORMULATION
  
% Model Parameters
u     =  90;                            % unemployment benefit
v     =  95;                            % benefit of pure leisure
wbar  = 100;                           	% long-run mean wage
gamma = 0.40;                           % wage reversion rate
p0    = 0.20;                         	% probability of finding job
p1    = 0.90;                           % probability of keeping job
sigma = 5;                             	% standard deviation of wage shock
delta = 0.95;                           % discount factor
   
% Continuous State Shock Distribution
m = 15;                                	% number of wage shocks
[e,w] = qnwnorm(m,0,sigma^2);           % wage shocks

% Stochastic Discrete State Transition Probabilities
q = zeros(2,2,2);
q(1,2,2) = p0;
q(2,2,2) = p1;
q(:,1,:) = 1-q(:,2,:);

% Model Structure
model.func = @func;                     % model functions
model.params = {u v wbar gamma};      	% function parameters
model.discount = delta;                	% discount factor
model.ds = 1;                           % dimension of continuous state
model.dx = 0;                           % dimension of continuous action
model.ni = 2;                           % number of discrete states
model.nj = 2;                           % number of discrete actions
model.e  = e;                          	% unit profit contribution shocks
model.w  = w;                          	% unit profit shock probabilities
model.q  = q;                          	% discrete state transition probabilities

% Approximation Structure
n = 150;                                % number of collocation nodes
wmin =   0;                             % minimum wage
wmax = 200;                             % maximum wage
basis = fundefn('spli',n,wmin,wmax);    % basis functions


%% SOLUTION

% Solve Bellman Equation
[c,sr,vr,xr,resid] = dpsolve(model,basis); 

% Compute and Print Critical Action Wages
wcrit1 = interp1(vr(:,1,1)-vr(:,1,2),sr,0);
vcrit1 = interp1(sr,vr(:,1,1),wcrit1);
fprintf('Critical Search Wage = %5.1f\n'    ,wcrit1) 
wcrit2 = interp1(vr(:,2,1)-vr(:,2,2),sr,0);
vcrit2 = interp1(sr,vr(:,2,1),wcrit2);
fprintf('Critical Quit Wage   = %5.1f\n\n'  ,wcrit2)

% Plot Action-Contingent Value Function - Unemployed
figure
hold on
plot(sr,squeeze(vr(:,1,:)))
yl = ylim;
plot([wcrit1 wcrit1],[yl(1) vcrit1],'k--')
bullet(wcrit1,yl(1))
bullet(wcrit1,vcrit1)
ftext(wcrit1,yl(1),'\bf w_0^*','left','bottom')
legend('Do Not Search','Search','Location','NW')
xlabel('Wage');
ylabel('Action-Contingent Value, Unemployed')

% Plot Action-Contingent Value Function - Employed
figure
hold on
plot(sr,squeeze(vr(:,2,:)))
yl = ylim;
plot([wcrit2 wcrit2],[yl(1) vcrit2],'k--')
bullet(wcrit2,yl(1))
bullet(wcrit2,vcrit2)
ftext(wcrit2,yl(1),'\bf w_1^*','left','bottom')
legend('Quit','Work','Location','NW')
xlabel('Wage');
ylabel('Action-Contingent Value, Employed')

% Plot Residual
figure
hold on
plot(sr,100*resid./max(vr,[],3))
plot(sr,0*sr,'k--')
legend('Idle','Active')
title('Bellman Equation Residual')
xlabel('Wage')
ylabel('Percent Residual')


%% SIMULATION

% Simulate Model
rand('seed',0.945);
nper = 41; 
nrep = 10000;
sinit = wbar*ones(nrep,1);
iinit = ones(nrep,1);
[ssim,xsim,isim] = dpsimul(model,basis,nper,sinit,iinit,sr,vr);

% Compute Ergodic Moments
isim = isim-1;
savg = mean(ssim(:,nper));
iavg = mean(isim(:,nper));
sstd = std(ssim(:,nper));
istd = std(isim(:,nper));

% Print Ergodic Moments
fprintf('Ergodic Means\n') 
fprintf('   Wage        = %5.1f\n'    ,savg)
fprintf('   Employment  = %5.2f\n\n'  ,iavg)
fprintf('Ergodic Standard Deviations\n') 
fprintf('   Wage        = %5.1f\n'    ,sstd)
fprintf('   Employment  = %5.2f\n\n'  ,istd)

% Plot Expected Discrete State Path
figure
hold on
plot(0:nper-1,mean(isim),'k')
title('Probability of Employment')
xlabel('Period')
ylabel('Probability')

% Plot Simulated and Expected Continuous State Path
figure
hold on
plot(0:nper-1,ssim(1:3,:))
plot(0:nper-1,mean(ssim),'k')
title('Simulated and Expected Wage')
xlabel('Period')
ylabel('Wage')


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

function out = func(flag,w,x,i,j,in,e,u,v,wbar,gamma)

switch flag
case 'f'
  out = (j==1)*v + (j==2)*((i==1)*u+(i==2)*w);
case 'g'
  out = wbar+gamma*(w-wbar)+e;
end
%%