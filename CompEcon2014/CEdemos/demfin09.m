function demfin09

%% DEMFIN09 Financial Asset Calibration

% Preliminary tasks
demosetup(mfilename)


%% FORMULATION

% Model Parameters
kappa  =  0.0363;
alpha  =  0.0692; 
sigma  =  0.0272;
 
% Model Structure
model.func     = @func;
model.T        = 30;
model.american = 0;
model.params   = {kappa,alpha,sigma};

% Approximation Structure
n = 20;
basis = fundefn('cheb',n,0,2);
s = funnode(basis);


%% SOLUTION

% Solve Model
optset('finsolve','keepall',1);
c = finsolve(model,basis,'lines',s,120);


%% ANALYSIS

% Calibrate to the data
y = [4.44 4.49 4.51 4.63 4.63 4.62 4.82 4.77 5.23];
tau = [.25 .5 1 2 3 5 7 10 30];
V = exp(-y/100.*tau);
t = (0:0.25:30);
tind = [2 3 5 9 13 21 29 41 121];
s = findstate(c(:,tind),basis,alpha,V);

% Create plots
Vhat = funeval(c,basis,s);
yhat = -100*log(Vhat)./t; yhat(1) = 100*s;

figure
plot(tau,y,'*',t,yhat)
title('Actual and Fitted Bond Yields')
xlabel('Time to Maturity')
ylabel('Yield')

disp('Model short interest rate and 3-month rate')
disp([s*100 y(1)])

% Save Plots as EPS Files
printfigures(mfilename,1)


%% FUNCTION FILE

function out = func(flag,S,kappa,alpha,sigma)
n = size(S,1);
switch flag
case 'rho'
  out = S;
case 'mu'
  out =  kappa*(alpha-S);
case 'sigma'
  out = sigma*sqrt(S);
case 'delta'
  out = [];
case 'V0'
  out = ones(n,1);
end