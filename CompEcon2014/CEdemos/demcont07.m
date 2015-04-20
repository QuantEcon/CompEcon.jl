function demcont07

%% DEMCONT07 Simple log-linear example in lecture notes

% Preliminary tasks
demosetup(mfilename)


%% FORMULATION

% Model Parameters
alpha =  0.2;       % capital share
rho   =  0.05;      % discount rate
sigma =  0.2 ;      % volatility

b = 1/rho;
a = log(rho)/rho + (alpha-rho-0.5*sigma^2)/(rho^2);


%% SOLVE BELLMAN EQUATION DIRECTLY USING BROYDEN

% Approximation Structure
n = 35;
smin = 0.2;
smax = 0.5;
basis = fundefn('cheb',n,smin,smax);
s = funnode(basis);

% Solve Bellman Equation Directly Using BROYDEN
v = a + b*log(s);
v = s;
c = funfitxy(basis,s,v);
optset('broyden','showiters',1)
c = broyden(@resid,c,basis,alpha,rho,sigma,s);

% Compute Residual
s = nodeunif(1000,smin,smax);
v = funeval(c,basis,s);
x = 1./funeval(c,basis,s,1);
res = resid(c,basis,alpha,rho,sigma,s);


% %% SOLVE BELLMAN EQUATION USING SCSOLVE
% 
% % Approximation Structure
% n = 15;
% smin = 0.2;
% smax = 0.3;
% basis = fundefn('cheb',n,smin,smax);
% s = funnode(basis);
% 
% % Model Structure
% model.func   = @func;
% model.rho    = rho;
% model.params = {alpha,sigma};
% 
% % Solve Bellman Equation Using SCSOLVE
% v = a + b*log(s);
% [c,s,v,x,res] = scsolve(model,basis,v);


%% PLOT SOLUTION

% Plot Value Function
figure
hold on
plot(s,v-a-b*log(s))
title('Value Function')
xlabel('Capital Stock')
ylabel('Value')

% Plot Shadow Price Function
figure
hold on
p = funeval(c,basis,s,1);
plot(s,p-b./s)
title('Shadow Price Function')
xlabel('Capital Stock')
ylabel('Shadow Price')

% Plot Optimal Policy
figure
hold on
plot(s,x-s/b)
title('Optimal Policy')
xlabel('Capital Stock')
ylabel('Consumption')

% Plot Residual
figure
plot(s,res,s,0*s,'k:')
title('Bellman Equation Residual')
xlabel('Capital Stock')
ylabel('Residual')


%% Save Plots as EPS Files
printfigures(mfilename,4)


%% Function File for scsolve
function out = func(flag,s,x,Vs,alpha,sigma)
switch flag
  case 'x'
    out = 1./Vs;
  case 'f'
    out = log(x);
  case 'g'
    out = alpha*s-x;
 case 'sigma'
   out = sigma*s;
end


%% Bellman Equation Residual
function res = resid(c,basis,alpha,rho,sigma,s)
V   = funeval(c,basis,s);
Vs  = funeval(c,basis,s,1);
Vss = funeval(c,basis,s,2);
x   = 1./Vs;
res = log(x) + Vs.*(alpha*s-x) + 0.5*Vss.*(sigma.*s).^2 - rho*V;