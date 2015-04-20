function demfin06

%% DEMFIN06 Compound Bond Option Pricing

% Preliminary tasks
demosetup(mfilename)


%% FORMULATION

% Model Parameters 
kappa  =  0.1;       % speed of mean reversion
alpha  =  0.05;      % long-run mean interest rate
sigma  =  0.1;       % interest rate volatility 
TB     =  30;        % bond maturity
K      =  0.2;       % exercise price
put    =  0;         % put indicator
TO     =  1;         % option maturity

% Model Structure
modelB.func     = @func;
modelB.T        = TB;
modelB.american = 0;
modelB.params   = {kappa,alpha,sigma};

% Approximation Structure
n = 20;
basisB = fundefn('cheb',n,0,2);


%% SOLUTION

% Solve Model
cB = finsolve(modelB,basisB);

% Create model variable the option
modelO.func     = @func;
modelO.T        = TO;
modelO.american = 0;
modelO.params   = {kappa,alpha,sigma,K,put,cB,basisB};

% Approximation Structure
n = 80;
basisO = fundefn('spli',n,0,2);

% Solve Model
cO = finsolve(modelO,basisO);


%% ANALYSIS

% Create plots
x = linspace(0,.25,201)';
Vhat = funeval(cO,basisO,x);

figure
plot(x,Vhat)
title('As a Function of Short Rate')
xlabel('r')
ylabel('Option Premium')
ylim([0 0.25])

BondVal = funeval(cB,basisB,x);
figure
plot(BondVal,Vhat)
title('As a Function of Bond Value')
xlabel('Bond Value')
ylabel('Option Premium')
ylim([0 0.25])
xlim([0 .45])

% Save Plots as EPS Files
printfigures(mfilename,2)

%% FUNCTION FILE

function out = func(flag,S,kappa,alpha,sigma,K,put,cB,basisB)
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
  if nargin<6
    out = ones(n,1);
  else
    bondval = funeval(cB,basisB,S);
    if put
      out = max(K-bondval,0);
    else
      out = max(bondval-K,0);
    end  
  end
end