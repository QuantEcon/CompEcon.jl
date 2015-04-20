function demfin01

%% DEMFIN01 Cox-Ingersoll-Ross Bond Pricing Model

% Preliminary tasks
demosetup(mfilename)


%% FORMULATION

% Model Parameters
T = 30;
kappa = .1;
alpha = .05;
sigma = 0.1;

% Model Structure
model.func     = @func;
model.T        = T;
model.american = 0;
model.params   = {kappa,alpha,sigma};

% Approximation Structure
n = 20;
basis = fundefn('cheb',n,0,2);


%% SOLUTION

% Solve Model
s = funnode(basis);
c = finsolve(model,basis,'lines',s,1);

% Compute exact solution and plot solution and errors
sigma2 = sigma*sigma;
gam = sqrt(kappa.^2+2*sigma2);
egam = exp(gam*T)-1;
denom = (gam+kappa).*egam+2*gam;
B = 2*egam./denom;
A = (2*gam*exp((kappa+gam)*T/2)./denom).^(2*kappa*alpha./sigma2);


%% ANALYSIS

x = linspace(0,0.25,101)';
V = A*exp(-B*x);
Vhat = funeval(c(:,end),basis,x);

figure
plot(x,Vhat)
title('Bond Price')
xlabel('r')
ylabel('V(r)')

figure
plot(x,V-Vhat)
title('Approximation Errors')
xlabel('r')
ylabel('Error')

fprintf('Maximum error on [0,0.25]: %10.4e\n',max(abs(V-Vhat)))

% Save Plots as EPS Files
printfigures(mfilename,2)


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