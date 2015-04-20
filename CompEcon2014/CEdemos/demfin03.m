function demfin03

%% DEMFIN03 Heston's stochastic volatility option pricing model

% Preliminary tasks
demosetup(mfilename)


%% FORMULATION

% Model Parameters
r      =  0.05;     % risk free interest rate
delta  =  0;        % rate of dividend payments
kappa  =  1;        % mean-reversion parameter on volatility 
m      =  0.05;     % long-run mean volatility
sigma  =  0.2;      % volatility of volatility
rho    =  -0.5;     % correlation between price and volatility
K      =  1;        % strike price
put    =  1;        % 0 = call, 1 = put
T      =  1;        % time to maturity

% Model Structure
model.func     = @func;
model.T        = T;
model.american = 0;
model.params   = {r,delta,kappa,m,sigma,rho,K,put};

% Approximation Structure
n = [100 10];
smin = [log(0.01*K) 0.1*m];
smax = [log(5*K) 4*m];
basis = fundefn('spli',n,smin,smax);
s = funnode(basis);


%% SOLUTION

% Solve Model
N = 1000;
c = finsolve(model,basis,'implicit',s,N);

% Compute BS solution and produce plots
p = linspace(log(.5*K),log(1.5*K),251)';
nu = m;
S = exp(p);

V1 = funeval(c(:,end),basis,{p,nu});
V2 = bs(sqrt(nu),S,K,r,delta,T,put);


%% ANALYSIS

figure
plot(S,V1,S,V2)
title('Option Values')
xlabel('S')
ylabel('Premium')
xlim([.75*K,1.5*K])
legend('SV','Black-Scholes')

isigma = impvol(V1,S,K,r,delta,T,put,0);
figure
plot(1./S,isigma,1./S,sqrt(nu)+zeros(size(S,1),1))
title('Implied Volatilities')
xlabel('K')
ylabel('Volatility')
xlim([0.75*K,1.5*K])
legend('SV','Black-Scholes')

nu = [.005 .05 .1 .125 .15 .175 .2]';
S = gridmake(p,nu);
V = reshape(funeval(c(:,end),basis,S),251,7);
figure
plot(exp(p),V)
title('Option Values for Alternative Values of \nu')
xlabel('S')
ylabel('Premium')
nn = length(nu);
legend([repmat('\nu  =  ',nn,1) reshape(sprintf('%5.3f',nu),5,nn)'],1)

% Save Plots as EPS Files
printfigures(mfilename,3)


%% FUNCTION FILE

function out = func(flag,S,r,delta,kappa,m,sigma,rho,K,put)
n = size(S,1);
switch flag
case 'rho'
  out = r+zeros(n,1);
case 'mu'
  out =  [r-delta-0.5*S(:,2)  kappa*(m-S(:,2))];
case 'sigma'
  out = sqrt(S(:,2));
  out = [out rho*sigma*out zeros(n,1) sqrt(1-rho*rho)*sigma*out];
case 'delta'
  out = [];
case 'V0'
  if put
    out = max(0,K-exp(S(:,1)));
  else
    out = max(0,exp(S(:,1))-K);
  end
end