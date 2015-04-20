function demfin05

%% DEMFIN05 Barrier Option Pricing

% Preliminary tasks
demosetup(mfilename)


%% FORMULATION

% Moldel Parameters
r       =  0.05;               % risk free interest rate
delta   =  0;                  % dividend rate 
sigma   =  0.2;                % volatility  
K       =  1;                  % exercise price 
put     =  0;                  % put indicator
T       =  1;                  % time-to-maturity
Sb      =  0.8;                % barrier 

% Model Structure
model.func     = @func;
model.T        = T;
model.american = 0;
model.params   = {r,delta,sigma,K,put};
model.barrier  = [0 Sb 1;0 inf 1];

% Approximation Structure
n = 75; 
basis = fundefn('spli',n,Sb,2*K);    


%% SOLUTION

% Solve Model
c = finsolve(model,basis);


%% ANALYSIS

% Create plots
S = sort([linspace(.5*K,2*K,501)';Sb-eps;Sb]);
Vhat = funeval(c,basis,S);
Vhat(S<Sb) = 0;
Vbs = bs(sigma,S,K,r,delta,T,put);
V = Vbs-(S./Sb).^(1-2*r/sigma^2).*bs(sigma,Sb^2./S,K,r,delta,T,put);
V(S<Sb) = 0;

figure
plot(S,Vhat,S,Vbs)
title('Down-and-out Call Option Premium')
xlabel('S')
ylabel('Pemium')
legend({'Barrier Option','Vanilla Option'},2)
xlim([.5 1.5]*K)
ylim([0 0.5])

figure
plot(S,V-Vhat)
title('Approximation Error')
xlabel('S')
ylabel('Error')
xlim([Sb,2*K])

disp('Sb Vb(Sb) Vbs(Sb)')
disp([S(S == Sb) Vhat(S == Sb) Vbs(S == Sb)])

n = 75; 
lambda = .5;
Delta = (2*K-Sb)/(n-2);
basis = fundefn('lin',n,Sb-lambda*Delta,2*K);    
% Call solution algorithm
s = funnode(basis);
s(abs(s-Sb)<1e-13) = Sb;
c = finsolve(model,basis,[],s);

Vhat = funeval(c,basis,S);
figure
plot(S,V-Vhat)
title('Approximation Error with Barrier Not a Node')
xlabel('S')
ylabel('Error')
xlim([Sb,2*K])

% Save Plots as EPS Files
printfigures(mfilename,3)


%% FUNCTION FILE

function out = func(flag,S,r,delta,sigma,K,put)
switch flag
case 'rho'
  n = size(S,1);
  out = r+zeros(n,1);
case 'mu'
  out =  (r-delta)*S;
case 'sigma'
  out = sigma*S;
case 'delta'
  out = [];
case 'V0'
  if put
    out = max(0,K-S);
  else
    out = max(0,S-K);
  end
end