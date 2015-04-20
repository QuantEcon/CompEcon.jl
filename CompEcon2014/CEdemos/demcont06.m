function demcont06

%% DEMCONT06 Odd Nonrenewable Resource Model
%
% Old demsc06

% Preliminary tasks
demosetup(mfilename)


%% FORMULATION

% Model Parameters
alpha = 0.15;
rho   = 0.05;


%% SOLVE USING SCSOLVE

% Approximation Structure
n = 3;
smin  = 0.01;
smax  = 100;
basis = fundefn('cheb',n,smin,smax);
s = funnode(basis);

% Model Structure
model.func   = @func;
model.rho    = rho;
model.params = {alpha};

% Solve Bellman Equation
[cv,s,v,x,resid] = scsolve(model,basis,s);


%% PLOT SOLUTION

% Plot Optimal Policy
figure
plot(s,x)
title('Optimal Policy')
xlabel('Resource Stock')
ylabel('Harvest')

% Plot Value Function
figure
plot(s,v)
title('Value Function')
xlabel('Resource Stock')
ylabel('Value')

% Plot Shadow Price Function
figure
plot(s,funeval(cv,basis,s,1))
title('Shadow Price Function')
xlabel('Resource Stock')
xlabel('Shadow Price')

% Plot Residual
figure
plot(s,resid,s,0*resid,'k:')
title('Bellman Equation Residual')
xlabel('Resource Stock')
ylabel('Residual')


%% Save Plots as EPS Files
printfigures(mfilename,4)


%% Function File for scsolve
function out = func(flag,s,x,Vs,alpha)
switch flag
  case 'x'
    out = (s.^(-alpha/(1-alpha)).*Vs).^(-1/alpha);
  case 'f'
    out = x.^(1-alpha);
  case 'g'
    out = -x.*(1-alpha).*(s.^(-alpha/(1-alpha)));
  case 'sigma'
    out = [];
end