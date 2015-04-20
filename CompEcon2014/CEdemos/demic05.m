function demic05

%% DEMIC05 Cash Management Model

% Preliminary tasks
demosetup(mfilename)


%% Formulation

% Define parameters
mu    = 0;
sigma = 0.2;
rho   = 0.04;
r     = 0.05;
c     = 0.75;
C     = 0.25;
f     = 1;
F     = 0.5;

% Create model variable
model.func=@func;
model.params={mu,sigma,rho,r,c,C};
model.xindex=[1 1;2 1];
model.F=[f;F];

% Define starting values
x=[0 1;2 1];

% Call solver
n=10;
optset('broyden','showiters',1)
[cv,basis,x]=icsolve(model,x,n);

% Plot results
S=linspace(0,x(2,1),101)';
V=funeval(cv,basis,S);
Vstar=funeval(cv,basis,x(:));
figure
plot(S,V,'k',x(:),Vstar,'k*')
title('Value Function')
xlabel('S')
ylabel('V')

dV=funeval(cv,basis,S,1);
dVstar=funeval(cv,basis,x(:),1);
figure
plot(S,dV,'k',[0;3],ones(2,1)*(r/rho+[c -C]),'k-',x(:),dVstar,'k*')
title('Marginal Value Function')
xlabel('S')
ylabel('V''')

% Display selected values
format short
disp('         x      V(x)     V''(x)')
disp([x([1;3;4;2]) Vstar([1;3;4;2]) dVstar([1;3;4;2])])

% Save Plots as EPS Files
printfigures(mfilename,2)


% Model function file for cash management example
function out=func(flag,s,mu,sigma,rho,r,c,C)
switch flag
  case 'f'
    out=zeros(size(s,1),1);
  case 'mu'
    out=mu+zeros(size(s,1),1);
  case 'sigma'
    out=sigma+zeros(size(s,1),1);
  case 'rho'
    out=rho+zeros(size(s,1),1);
  case 'R+'
    out=[(-c-r/rho)*(s(2)-s(1)) (c+r/rho) (-c-r/rho) 0];
  case 'R-'
    out=[(r/rho-C)*(s(1)-s(2)) (r/rho-C) (C-r/rho) 0];
end