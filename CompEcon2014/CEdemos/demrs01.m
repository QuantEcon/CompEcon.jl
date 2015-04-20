function demrs01

%% DEMRS01 Asset Abandonment Model

% Preliminary tasks
demosetup(mfilename)


%% Formulation

% Define parameters
c     = 0.5;
mu    = 0;
sigma = 0.2;
rho   = 0.1;

% Pack model structure
model.func=@func;
model.params={c,mu,sigma,rho};
model.xindex=[1 0 1 2 0;1 0 0 0 1];

% Set approximation size and initial values
n=101;
x=[0;20];

% Call solver
[cv,basis,x]=rssolve(model,x,n,'cheb');
cv=cv{1};
basis=basis{1};

% Plot results
S=linspace(0,1,1001)';
V=funeval(cv,basis,S);

sstar=x(1);

V(S<sstar)=0;

figure
plot(S,V,'k',sstar,funeval(cv,basis,sstar),'k*')
title('Value Function')
xlabel('P')
ylabel('V');

figure
Vs=funeval(cv,basis,S,1);
Vs(S<sstar)=0;
plot(S,Vs,'k',sstar,funeval(cv,basis,sstar,1),'k*')
title('Marginal Value Function')
xlabel('P')
ylabel('V''')

S=linspace(0,basis.b,1001)';
V=funeval(cv,basis,S);
Vs=funeval(cv,basis,S,1);
Vss=funeval(cv,basis,S,2);

sstar=x(1);

V(S<sstar)=0;
Vs(S<sstar)=0;
Vss(S<sstar)=0;

% Exact solution
beta=roots([sigma^2/2 mu-sigma^2/2 -rho]);
beta=beta(beta<0);
pstar=(rho-mu)*beta*c/rho/(beta-1);
A=-pstar.^(1-beta)/(rho-mu)/beta;
VV=S/(rho-mu)-c/rho+A*S.^beta;
VV(S<pstar)=0;

e=rho*V-(S-c).*(S>sstar)-mu*S.*Vs-0.5*sigma^2*S.^2.*Vss;

figure
plot(S,e)
title('Aproximation Residual')
xlabel('P')
ylabel('Residual')

figure
plot(S,VV-V(:,end))
title('Approximation Error')
xlabel('P')
ylabel('Error')

disp('P*: Approx    Exact     Error')
disp([sstar pstar sstar-pstar])

% Save Plots as EPS Files
printfigures(mfilename,4)


function out=func(flag,s,x,c,mu,sigma,rho)
switch flag
  case 'f'
    out=(s-c);
  case {'g','mu'}
    out=mu*s;
  case 'sigma'
    out=sigma*s;
  case 'rho'
    out=rho+zeros(size(s,1),1);
  case 'reward'
    out=[0 0 0;s(2)/(rho-mu)-c/rho 1/(rho-mu) 0];
  otherwise
    error('invalid flag')
end