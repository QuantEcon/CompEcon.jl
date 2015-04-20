function demic03

%% DEMIC03 Storage Management Model

% Preliminary tasks
demosetup(mfilename)


%% Formulation

% Define parameters
mu    = -0.2;
sigma = 0.05;
k     = 0.05;
P     = 2;
F     = 3;
rho   = 0.1;

% Create model variable
model.func=@func;
model.params={mu,sigma,k,P,rho};
model.xindex=[1 1;0 0];
model.F=[F;0];

% Define starting values
smax=8;
x=[0 smax;smax smax];

% Call solver
n=10;
optset('broyden','showiters',1)
[cv,basis,x]=icsolve(model,x,n);
Sstar=x(1,2);

% Plot results
S=linspace(0,smax,101)';
V=funeval(cv,basis,S);
Vstar=funeval(cv,basis,Sstar);
figure
plot(S,V,'k',Sstar,Vstar,'k*',S,Vstar+P*(S-Sstar))
title('Value Function')
xlabel('S')
ylabel('V')
text(12,19,['S^* =' num2str(Sstar)])

disp('     S         V')
disp([0 V(1);Sstar Vstar])

% Save Plots as EPS Files
printfigures(mfilename,1)


% Model function file for storage example
function out=func(flag,s,mu,sigma,k,P,rho)
switch flag
  case 'f'
    out=-k*s;
  case 'mu'
    out=mu+zeros(size(s,1),1);
  case 'sigma'
    out=sigma+zeros(size(s,1),1);
  case 'rho'
    out=rho+zeros(size(s,1),1);
  case 'R+'
    out=[P*(s(1)-s(2)) P -P 0];
  case 'R-'
    out=zeros(1,4);
end