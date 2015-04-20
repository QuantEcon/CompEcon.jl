function demdp01

%% DEMDP01 Timber Harvesting Model - Cubic Spline Approximation
%
% Profit maximizing owner of a commercial tree stand must decide when to
% clearcut the stand.
%
% States
%     s       stand biomass
% Actions
%     j       clear cut (2) or do not clear cut (1)
% Parameters
%     price   unit price of biomass
%     kappa   clearcut-replant cost
%     smax    stand carrying capacity
%     gamma   biomass growth parameter
%     delta   discount factor

% Preliminary tasks
demosetup(mfilename)


%% FORMULATION

% Model Parameters
price = 1.0;                           	% unit price of biomass
kappa = 0.2;                            % clearcut-replant cost
smax  = 0.5;                            % stand carrying capacity
gamma = 0.1;                            % biomass growth parameter
delta = 0.9;                            % discount factor

% Approximation Structure
n = 200;                                % number of collocation nodes
basis = fundefn('spli',n,0,smax);       % basis functions
  
% Model Structure
model.func = @func;                     % model functions
model.params = {price kappa smax gamma};% function parameters
model.discount = delta;                	% discount factor
model.ds = 1;                           % dimension of continuous state
model.dx = 0;                           % dimension of continuous action
model.ni = 0;                           % number of discrete states
model.nj = 2;                           % number of discrete actions


%% SOLUTION

% Solve Bellman Equation
[c,sr,vr,xr,resid] = dpsolve(model,basis); 

% Plot Action-Contingent Value Functions
figure
hold on
plot(sr,vr)
legend('Grow','Clear-Cut','Location','N')
title('Action-Contingent Value Functions')
xlabel('Biomass')
ylabel('Value of Stand')

% Compute and Plot Optimal Harvesting Stock Level
scrit = interp1(vr(:,1)-vr(:,2),sr,0);
vcrit = interp1(sr,vr(:,1),scrit);
plot([scrit scrit],[-0.2 vcrit],'k--')
bullet(scrit,-0.2)
bullet(scrit,vcrit)
ftext(scrit,-0.2,'\bf s^*','left','bottom')
fprintf('Optimal Biomass Harvesting Level = %5.2f\n',scrit)

% Plot Residual
figure
hold on
plot(sr,100*resid./max(vr,[],2))
plot(sr,0*sr,'k--')
title('Bellman Equation Residual')
xlabel('Biomass')
ylabel('Percent Residual')


%% SIMULATION

% Simulate Model
nper = 31;      % Number of periods simulated
time = 0:nper-1;% Periods simulated
sinit = 0;      % Initial value of continuous state
iinit = 1;      % Initial value of discrete state
[ssim,xsim,isim,jsim] = dpsimul(model,basis,nper,sinit,iinit,sr,vr);

% Compute Optimal Rotation Cycle
fprintf('Optimal Rotation Cycle           = %5i\n',min(time(jsim==2)))

% Plot State Path
figure
plot(time,ssim)
xlabel('Period')
ylabel('Biomass')


%% Save Plots as EPS Files
printfigures(mfilename,3)


%% DPSOLVE FUNCTION FILE
%
%    User-supplied function called by dpsolve that returns the bound,
%    reward, and continuous state transition function values and
%    derivatives with respect to the continuous action x at an arbitrary
%    number ns of states and actions aaccording to the format
%      [out1,out2,out3] = func(flag,s,x,i,j,in,e,<params>)
%    For further information regarding input format, type > help dpsolve or
%    reference the DPSOLVE Handbook.

function out = func(flag,s,x,i,j,in,e,price,kappa,smax,gamma)
n = length(s);
switch flag
  case 'f'
    out = (price*s-kappa).*(j-1);
  case 'g'
    if j==1
      out = s+gamma*(smax-s);
    else
      out = gamma*smax*ones(n,1);
    end
end
%%