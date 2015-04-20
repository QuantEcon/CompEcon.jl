function demdp10

%% DEMDP10 Water Resource Management Model
%
% Public authority must decide how much water to release from a reservoir so 
% as to maximize benefits derived from agricultural and recreational uses.
%
% States
%     s       reservoiur level at beginning of summer
% Actions
%     x       quantity of water released for irrigation
% Parameters
%     a       producer benefit function parameters
%     b       recreational user benefit function parameters
%     delta   discount factor

% Preliminary tasks
demosetup(mfilename)


%% FORMULATION

% Model Parameters
a = [1 -2];                                 % producer benefit function parameters
b = [2 -3];                                 % recreational user benefit function parameters
ymean = 1.0;                                % mean rainfall
sigma = 0.2;                                % rainfall volatility
delta = 0.9;                                % discount factor

% Continuous State Shock Distribution
m = 3;                                      % number of rainfall shocks
[e,w] = qnwlogn(m,log(ymean)-sigma^2/2,sigma^2); % rainfall shocks and proabilities

% Model Structure
model.func = @func;                         % model functions
model.params = {a b};                       % function parameters
model.discount = delta;                     % discount factor
model.ds = 1;                               % dimension of continuous state
model.dx = 1;                               % dimension of continuous action
model.ni = 0;                               % number of discrete states
model.nj = 0;                               % number of discrete actions
model.e  = e;                               % continuous state shocks
model.w  = w;                               % continuous state shock probabilities

% Approximation Structure
n    = 15;                                  % number of collocation nodes
smin =  2;                                  % minimum state
smax =  8;                                  % maximum state
basis = fundefn('cheb',n,smin,smax);        % basis functions
s = funnode(basis);                         % collocaton nodes
  
% Deterministic Steady-State
xstar = 1;                                  % deterministic steady-state action
sstar = 1+(a(1)*(1-delta)/b(1))^(1/b(2));   % deterministic steady-state stock

% Check Model Derivatives
dpcheck(model,sstar,xstar);


%% SOLUTION

% Compute Linear-Quadratic Approximation at Collocation Nodes
[vlq,xlq] = lqapprox(model,s,sstar,xstar); 

% Solve Bellman Equation
[c,s,v,x,resid] = dpsolve(model,basis,vlq,xlq);

% Plot Optimal Policy
figure
plot(s,x)
title('Optimal Irrigation Policy')
xlabel('Reservoir Level')
ylabel('Irrigation')

% Plot Value Function
figure
plot(s,v)
title('Value Function')
xlabel('Reservoir Level')
ylabel('Value')

% Plot Shadow Price Function
figure
pr = funeval(c,basis,s,1);
plot(s,pr)
title('Shadow Price Function')
xlabel('Reservoir Level')
ylabel('Shadow Price')

% Plot Residual
figure
hold on
plot(s,resid)
plot(s,0*s,'k--')
title('Bellman Equation Residual')
xlabel('Reservoir Level')
ylabel('Residual')


%% SIMULATION

% Simulate Model
rand('seed',0.945);
nper = 31;
nrep = 100000;
sinit = smin*ones(nrep,1);
iinit = ones(nrep,1);
[ssim,xsim] = dpsimul(model,basis,nper,sinit,iinit,s,v,x);

% Plot Simulated State Path
figure
hold on
plot(0:nper-1,ssim(1:3,:))
plot(0:nper-1,mean(ssim),'k')
title('Simulated and Expected Reservoir Level')
xlabel('Year')
ylabel('Reservoir Level')

% Plot Simulated Policy Path
figure
hold on
plot(0:nper-1,xsim(1:3,:))
plot(0:nper-1,mean(xsim),'k')
title('Simulated and Expected Irrigation')
xlabel('Year')
ylabel('Irrigation')

% Compute Ergodic Moments
ssim = ssim(:,nper);
xsim = xsim(:,nper);
savg = mean(ssim(:)); 
xavg = mean(xsim(:)); 
sstd = std(ssim(:)); 
xstd = std(xsim(:)); 

% Print Steady-State and Ergodic Moments
fprintf('Deterministic Steady-State\n') 
fprintf('   Reservoir Level    = %5.2f\n'    ,sstar)
fprintf('   Irrigation         = %5.2f\n\n'  ,xstar)
fprintf('Ergodic Means\n') 
fprintf('   Reservoir Level    = %5.2f\n'    ,savg)
fprintf('   Irrigation         = %5.2f\n\n'  ,xavg)
fprintf('Ergodic Standard Deviations\n') 
fprintf('   Reservoir Level    = %5.2f\n'    ,sstd)
fprintf('   Irrigation         = %5.2f\n\n'  ,xstd)

% Compute and Plot Ergodic State Distribution
[qq,ss] = ksdensity(ssim(:),'support','positive');
figure
plot(ss,qq)
title('Ergodic Reservoir Level Distribution')
xlabel('Reservoir Level')
ylabel('Probability')
xlim([2 6])  


%% Save Plots as EPS Files
printfigures(mfilename,7)


%% DPSOLVE FUNCTION FILE
%
%    User-supplied function called by dpsolve that returns the bound,
%    reward, and continuous state transition function values and
%    derivatives with respect to the continuous action x at an arbitrary
%    number ns of states and actions aaccording to the format
%      [out1,out2,out3] = func(flag,s,x,i,j,in,e,<params>)
%    For further information regarding input format, type > help dpsolve or
%    reference the DPSOLVE Handbook.

function [out1,out2,out3] = func(flag,s,x,i,j,in,e,a,b)

switch flag
case 'b'  
  out1 = zeros(size(s));                 
  out2 = s;        
  out3 = [];                      
case 'f'
  out1 = (a(1)/(1+a(2)))*x.^(1+a(2))+(b(1)/(1+b(2)))*(s-x).^(1+b(2));   
  out2 = a(1)*x.^a(2)-b(1)*(s-x).^b(2);                                 
  out3 = a(1)*a(2)*x.^(a(2)-1)+b(1)*b(2)*(s-x).^(b(2)-1);               
case 'g' 
  out1 = s-x+e;                          
  out2 = -ones(size(s));                 
  out3 = zeros(size(s)); 
end 
%%