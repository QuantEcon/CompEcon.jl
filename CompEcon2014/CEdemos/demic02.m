function demic02

%% DEMIC02 Timber Harvesting Model

% Preliminary tasks
demosetup(mfilename)


%% Formulation

% Define parameters
alpha = 0.1;
m     = 1;
sigma = 0.05;
P     = 3;
C     = 0.15;
rho   = 0.1;

% Create model variable
model.func=@func;
model.params={alpha,m,sigma,P,rho};
model.xindex=[0 0;2 0];
model.F=[0;C];

% Define starting values
x=[0 0;0.5 0];

% Call solver
n=15;
optset('broyden','showiters',1)
[cv,basis,x]=icsolve(model,x,n);
Sstar=x(2,1);

% Plot results
S=linspace(0,m/2,101)';
V=funeval(cv,basis,S);
Vstar=funeval(cv,basis,Sstar);
dV=funeval(cv,basis,S,1);
dVstar=funeval(cv,basis,Sstar,1);
V(S>Sstar)=Vstar+P*(S(S>Sstar)-Sstar);
dV(S>Sstar)=P;
figure
plot(S,V,'k',Sstar,Vstar,'k*')
title('Value Function')
xlabel('S')
ylabel('V')
ylim([1.8 3.2])

figure
plot(S,dV,'k',Sstar,dVstar,'k*')
title('Marginal Value')
xlabel('S')
ylabel('V''')
ylim([1.8 3.2])

% Save Plots as EPS Files
printfigures(mfilename,2)


% Model function file for timber harvesting example
function out=func(flag,s,alpha,m,sigma,P,rho)
switch flag
  case 'f'
    out=zeros(size(s,1),1);
  case 'mu'
    out=alpha*(m-s);
  case 'sigma'
    out=sigma*sqrt(s);
  case 'rho'
    out=rho+zeros(size(s,1),1);
  case 'R+'
    out=zeros(1,4);
  case 'R-'
    out=[P*(s(1)-s(2)) P -P 0];
end