function demfin02

%% DEMFIN02 Black-Scholes option pricing model

% Preliminary tasks
demosetup(mfilename)


%% FORMULATION

% Model Parameters
r       =  0.05;               % risk free interest rate
deltaS  =  0;                  % dividend rate 
sigma   =  0.2;                % volatility  
K       =  1;                  % exercise price 
put     =  0;                  % put indicator
T       =  1;                  % time-to-maturity

% Model Structure
model.func     = @func;
model.T        = T;
model.params   = {r,deltaS,sigma,K,put};
model.american = 0;

% Approximation Structure
n = 51; 
basis = fundefn('lin',n,0,2*K);    
s = funnode(basis);


%% SOLUTION

% Solve Model
c = finsolve(model,basis,'lines',s);

% Compute exact solution and plot solution and errors
S = linspace(0,2*K,501)';
premium = bs(sigma,S,K,r,deltaS,T,put);


%% ANALYSIS

figure
plot(S,funeval(c,basis,S))
title('Call Option Premium')
xlabel('S')
ylabel('Premium')

figure
plot(S,premium-funeval(c,basis,S))
title('Approximation Errors')
xlabel('S')
ylabel('Error')

% Compute performance comparisons

stype = {'lines' 'explicit' 'implicit' 'CN' 'stiff'};
N = [1 75 75 75 75];
m = length(stype);
C = zeros(n,m);
tl = zeros(1,m);

basis = fundefn('lin',n,0,2*K);    
s = funnode(basis);

for i = 1:m
  tic;
  c = finsolve(model,basis,stype{i},s,N(i));
  tl(i) = toc;
  C(:,i) = c(:,end);
end
VL = funeval(C,basis,S);

ts = zeros(1,m);
N = [1 250 75 75 75];

basis = fundefn('spli',n,0,2*K);    
s = funnode(basis);
for i = 1:m
  tic;
  c = finsolve(model,basis,stype{i},s,N(i));
  ts(i) = toc;
  C(:,i) = c(:,end);
end
VS = funeval(C,basis,S);

disp('Maximum errors on [0,2K] for lines, explicit, implicit, CN and stiff methods')
disp('Piecewise Linear')
show(max(abs(premium(:,ones(5,1))-VL)),6)
disp('Cubic Spline')
show(max(abs(premium(:,ones(5,1))-VS)),6)

disp('Timing for lines, explicit, implicit, CN and stiff methods')
disp('Piecewise Linear')
show(tl)
disp('Cubic Spline')
show(ts)

% Save Plots as EPS Files
printfigures(mfilename,2)


%% FUNCTION FILE

function out = func(flag,S,r,deltaS,sigma,K,put)
switch flag
case 'rho' 
  out  =  r+zeros(size(S,1),1);
case 'mu'
  out  =  (r-deltaS)*S;
case 'sigma'
  out  =  sigma*S;
case 'delta'
 out  =  [];
case 'V0'
  if put
    out  =  max(0,K-S);
  else
    out  =  max(0,S-K);
  end
end