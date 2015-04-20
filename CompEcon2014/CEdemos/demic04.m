function demic04

%% DEMIC04 Capacity Choice Model

% Preliminary tasks
demosetup(mfilename)


%% Formulation

% Define parameters
P     = 2;
C     = 1;
delta = 0.5;
rho   = 0.1;

% Create model variable
model.func=@func;
model.params={P,C,delta,rho};
model.xindex=[2 0;0 0];
model.F=[0;0];

% Define starting values
xmax=20;
x=[0 0;xmax 0];

% Call solver
n=25;
optset('broyden','showiters',1)
[cv,basis,x]=icsolve(model,x,n);
Kstar=x(1,1);

% Plot results
K=linspace(0,xmax,101)';
V=funeval(cv,basis,K);
Vstar=funeval(cv,basis,Kstar);
V(K<Kstar)=Vstar-(Kstar-K(K<Kstar))*C;
figure
plot(K,V,'k',Kstar,Vstar,'k*')
title('Value Function')
xlabel('K')
ylabel('V')
text(3,Vstar,['K^* =' num2str(Kstar)])

disp('       K*     V(K*)')
disp([Kstar Vstar])

% Save Plots as EPS Files
printfigures(mfilename,1)


% Model function file for capacity choice example
function out=func(flag,s,P,C,delta,rho)
switch flag
  case 'f'
    out=P*log(s+1);
  case 'mu'
    out=-delta*s;
  case 'sigma'
    out=[];
  case 'rho'
    out=rho+zeros(size(s,1),1);
  case 'R+'
    out=[C*(s(1)-s(2)) C -C 0];
  case 'R-'
    out=zeros(1,4);
end