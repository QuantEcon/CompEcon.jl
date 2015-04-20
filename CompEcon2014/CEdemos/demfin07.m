function demfin07

%% DEMFIN07 Asian Option Pricing

% Preliminary tasks
demosetup(mfilename)


%% FORMULATION

% Model Parameters
r      =  0.1;
delta  =  0;
sigma  =  0.4;
L      =  1/2;
put    =  0;
tau    =  1/4;

% Model Structure
model.func     = @func;
model.T        = tau;
model.american = 0;
model.params   = {r,delta,sigma,L,put};

% Approximation Structure
n = 101; 
basis = fundefn('lin',n,0,1);


%% SOLUTION

% Solve Model
c = finsolve(model,basis);


%% ANALYSIS

% Transform solution to natural state space (y to S)
S = linspace(0.001,2,101)';
M = [.5 .75 1 1.25 1.5];
m = length(M);
Vhat = zeros(length(S),m);
for i = 1:m
  y = ((L-tau)*M(i))./S;
  Vhat(:,i) = S.*funeval(c,basis,y);
  if ~put, Vhat(y>basis.b,i) = 0; end
end

% Create plots
figure
plot(S,Vhat)
title('Call Option Premium')
xlabel('S')
ylabel(['V(S,M,' num2str(tau) ')'])
nn = length(M);
legend([repmat('M  =  ',nn,1) reshape(sprintf('%4.2f',M),4,nn)'],2)
set(gca,'ylim',max(0,get(gca,'ylim')))

y = linspace(basis.a,basis.b,101)';
figure
plot(y,funeval(c,basis,y))
title('Approximation Function')
xlabel('y')
ylabel('v(y,\tau)')
ylim([-0.1 0.6])

% Save Plots as EPS Files
printfigures(mfilename,2)


%% FUNCTION FILE

function out = func(flag,S,r,delta,sigma,L,put)
switch flag
case 'rho'
  n = size(S,1);
  out = delta+zeros(n,1);
case 'mu'
  out =  1-(r-delta)*S;
case 'sigma'
  out = sigma*S;
case 'delta'
  out = [];
case 'V0'
  if put
    out = max(0,S/L-1);
  else
    out = max(0,1-S/L);
  end
end