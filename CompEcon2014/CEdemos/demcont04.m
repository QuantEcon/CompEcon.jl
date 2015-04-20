function demcont04

%% DEMCONT04 Optimal Fish Harvest Model (Stochastic)
%
% Old demsc04

% Preliminary tasks
demosetup(mfilename)


%% FORMULATION

% Model Parameters
alpha = 0.5;
sigma = 0.5;
H     = 1;
P     = 1;
c     = 0.25;
rho   = 0.1;


%% SOLVE USING SCSOLVE

% Approximation Structure
n    = 1000;
smin = 0.01;
smax = 1.5;
basis = fundefn('spli',n,smin,smax);
s = funnode(basis);

% Model Structure
model.func   = @func;
model.rho    = rho;
model.params = {alpha,sigma,H,P,c};

% Solve Bellman Equation
[c,s,v,x,resid] = scsolve(model,basis,s);

% Approximate Optimal Switch Point and Value Functions
Vs = funeval(c,basis,s,1);
x  = feval(model.func,'x',s,[],Vs,model.params{:});
i  = sum(x==0);
sstar = (s(i)+s(i+1))/2;
S = [s(1:i);sstar;s(i+1:end)];
Vs = [Vs(1:i);funeval(c,basis,sstar,1);Vs(i+1:end)];


%% PLOT SOLUTION

% Plot Value Function
figure
hold on
plot(s,v)
vstar = funeval(c,basis,sstar);
plot(sstar,vstar,'k*','LineWidth',8)
title('Value Function')
xlabel('Fish Stock')
ylabel('Value')

% Plot Shadow Price Function
figure
hold on
p = funeval(c,basis,s,1);
plot(s,p)
pstar = funeval(c,basis,sstar,1);
plot(sstar,pstar,'k*','LineWidth',8)
xlim([smin smax])
title('Shadow Price Function')
xlabel('Fish Stock')
xlabel('Shadow Price')

% Plot Residual
figure
plot(s,resid,s,0*s,'k:')
title('Bellman Equation Residual')
xlabel('Fish Stock')
ylabel('Residual')


%% PLOT ERGODIC DISTRIBUTION

[cp,Ex] = itodensity(model,basis,c);
p = funeval(cp,basis,s);
basis0 = fundefn('cheb',51,basis.a,3);
s0 = linspace(basis0.a,basis0.b,201)';
cv0 = funfitxy(basis0,s0,(20*P)*s0);
[cp0,Ex0] = itodensity(model,basis0,cv0);
p0 = funeval(cp0,basis0,s0);

% Plot Ergodic Distribution
figure
plot(s,p,s0,p0)
title('Ergodic Distribution of Fish Stocks')
xlabel('Fish Stock')
ylabel('Probability')

disp(' ')
disp('      S*       E[S]    E[S|x=0]')
disp([sstar Ex Ex0])
disp(' ')
disp('Percentage of time inactive')
disp(funeval(cp,basis,sstar,-1))


%% Save Plots as EPS Files
printfigures(mfilename,4)


%% Function File for scsolve
function out = func(flag,s,x,Vs,alpha,sigma,H,P,c)
switch flag
  case 'x'
    out = H*(Vs<(P-c./s));
  case 'f'
    out = (P-c./s).*s.*x;
  case 'g'
    out = (alpha*(1-s)-x).*s;
  case 'sigma'
    out = sigma*s;
end