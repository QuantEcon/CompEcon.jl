function demdp00

%% DEMDP00 Timber Harvesting Model - Simple Linear Approximation
%
% Profit maximizing owner of a commercial tree stand must decide when to
% clear-cut the stand.  This program uses a simple linear approximation for
% the value function that can be derived by "hand:
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
price = 1.0;                           	% price of biomass
kappa = 0.2;                            % clearcut-replant cost
smax  = 0.5;                            % stand carrying capacity
gamma = 0.1;                            % biomass growth parameter
delta = 0.9;                            % discount factor


%% SOLUTION - Linear Value Function Approximant

% Solve collocation equation with two collocation nodes
s = [0.2;0.4]; 
c = zeros(2,1);
c = broyden(@residual,c,s,gamma,smax,delta,price,kappa);

% Compute residual and value function on refined grid
s = nodeunif(200,0,smax);
resid = residual(c,s,gamma,smax,delta,price,kappa);
v = [delta*(c(1)+c(2)*(s+gamma*(smax-s))) price*s-kappa+delta*(c(1)+c(2)*gamma*smax)];

% Plot Conditional Value Functions
figure
hold on
plot(s,v)
legend('Grow','Clear-Cut')
title('Conditional Value Functions')
xlabel('Biomass')
ylabel('Value of Stand')

% Compute and Plot Optimal Biomass Harvest Level
scrit = interp1(v(:,1)-v(:,2),s,0);
vcrit = interp1(s,v(:,1),scrit);
plot([scrit scrit],[-0.2 vcrit],'k--')
bullet(scrit,-0.2)
bullet(scrit,vcrit)
ftext(scrit,-0.2,'\bf s^*','left','bottom')
fprintf('Optimal Biomass Harvest Level = %5.2f\n',scrit) 

% Plot Value Function Residual
figure
hold on
plot(s,100*resid./max(v,[],2))
plot(s,0*s,'k--')
title('Bellman Equation Residual')
xlabel('Biomass')
ylabel('Percent Residual')


%% Save Plots as EPS Files
printfigures(mfilename,2)


%% Residual Function for n=2 Timber Harvesting Problem
function resid = residual(c,s,gamma,smax,delta,price,kappa)
h = s+gamma*(smax-s);
resid = c(1)+c(2)*s - max(delta*(c(1)+c(2)*h),price*s-kappa+delta*(c(1)+c(2)*gamma*smax));
%%        