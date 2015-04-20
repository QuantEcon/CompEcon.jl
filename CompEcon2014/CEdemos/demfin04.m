function demfin04

%% DEMFIN04 American Put Option Pricing

% Preliminary tasks
demosetup(mfilename)


%% FORMULATION

% Model Parameters
r       =  0.05;               % risk free interest rate
deltaS  =  0;                  % dividend rate 
sigma   =  0.2;                % volatility  
K       =  1;                  % exercise price 
put     =  1;                  % put indicator
T       =  1;                  % time-to-maturity

% Model Structure
model.func     = @func;
model.T        = T;
model.american = 1;
model.params   = {r,deltaS,sigma,K,put};

% Approximation Structure
n = 100; 
basis = fundefn('lin',n,0,2*K);    
s = funnode(basis);


%% SOLUTION

% Solve Model
N = 75;
optset('finsolve','keepall',1);
c = finsolve(model,basis,[],s,N);

% Compute Barone-Adesi/Whaley solution and plot solution and errors
S = linspace(basis.a,basis.b,301)';
V = funeval(c(:,end),basis,S);
[Vbaw,bs,sstarB] = baw(sigma,S,K,r,deltaS,T,put);


%% ANALYSIS

figure
plot(S,V,'k-',S,bs,'r-',sstarB,funeval(c(:,end),basis,sstarB),'k*')
title('Option Premium')
xlabel('S')
ylabel('Premium')
legend('American','European')

figure
plot(S,Vbaw-V,'k')
title('Approximation Error')
xlabel('S')
ylabel('Error')

% Compute and plot the optimal exercise boundary
V0 = feval(model.func,'V0',s,model.params{:});
temp = funeval(c,basis,s) == V0(:,ones(1,N+1)) & s(:,ones(1,N+1))<= K;
% use upper bound on S* at time values where the upper bound changes
sstar = s(sum(temp)'+1);
sstarl = s(sum(temp)');
tau = linspace(0,T,N+1)';
ind = [1;find(diff(sstar)~= 0)+1;N+1];
sstari = sstar(ind);
taui = tau(ind);
% end point adjustments
sstari(1) = K;
send = sstari(end-1)+(sstari(end-1)-sstari(end-2))/(taui(end-1)-taui(end-2))*(taui(end)-taui(end-1));
sstari(end) = (sstari(end)+send)/2;

figure
plot(taui,sstari,'k')
hold on; stairs(tau,sstar,'r--'); stairs(tau,sstarl,'r--'); hold off
title('Early Exercise Boundary')
xlabel('\tau')
ylabel('S^*')
ylim([0.8 1.05])

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