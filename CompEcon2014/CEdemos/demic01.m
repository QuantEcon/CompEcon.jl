function demic01

%% DEMIC01 Asset Replacement Model

% Preliminary tasks
demosetup(mfilename)


%% Formulation

% Define parameters
beta = [1; 0.05; -0.003];
P    = 2;
C    = 3;
rho  = 0.1;

% Create model variable
model.func=@func;
model.params={beta,P,rho};
model.xindex=[0 0;2 0];
model.F=[0;C];

% Define starting values
x=[0 0;100 0];

% Call solver
n=15;
optset('broyden','showiters',1)
[cv,basis,x]=icsolve(model,x,n);
Astar=x(2,1);

% Plot results
A=linspace(0,Astar,101)';
V=funeval(cv,basis,A);
figure
plot(A,V,'k',Astar,V(end),'k*')
title('Value Function')
xlabel('A')
ylabel('V')
h=text(15,19,['A^* =' num2str(Astar)]);
set(h,'HorizontalAlignment','right')
xlim([0 18])

% Save Plots as EPS Files
printfigures(mfilename,1)


% Model function file for asset replacement example
function out=func(flag,s,beta,P,rho)
switch flag
  case 'f'
    Q=(beta(3)*s+beta(2)).*s+beta(1);
    out=P*Q;
  case 'mu'
    out=ones(size(s,1),1);
  case 'sigma'
    out=[];
  case 'rho'
    out=rho+zeros(size(s,1),1);
  case 'R+'
    out=zeros(1,4);
  case 'R-'
    out=zeros(1,4);
end