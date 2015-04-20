function demrs03

%% DEMRS03 Dixit Entry/Exit Model
%
% ref.: Dixit & Pindyck, p. 218

% Preliminary tasks
demosetup(mfilename)


%% Formulation

% Define parameters
r     = 0.05;
mu    = 0;
sigma = 0.2;
C     = 1;
E     = 2;
I     = 5;

% Pack model structure
model.func=@func;
model.params={r,mu,sigma,C,E,I};
model.xindex=...
  [1 0 0 0 0;
  1 2 1 2 0;
  2 1 1 2 0;
  2 0 0 0 1];

% Set initial values and approximation size
x=[0; 3; 1; 50];
n=[15 85];

% Call solver
[cv,basis,x]=rssolve(model,x,n);

x=reshape(x,2,2);
pstar=[x(2,1);x(1,2)];

% Compute approximation value functions
P=linspace(0,3,501)';

V=[funeval(cv{1},basis{1},P),...
  funeval(cv{2},basis{2},P)];
V(P>x(2,1),1)=V(P>x(2,1),2)-I;
V(P<x(1,2),2)=V(P<x(1,2),1)-E;

dV=[funeval(cv{1},basis{1},P,1),...
  funeval(cv{2},basis{2},P,1)];
dV(P>x(2,1),1)=dV(P>x(2,1),2);
dV(P<x(1,2),2)=dV(P<x(1,2),1);
d2V=[funeval(cv{1},basis{1},P,2),...
  funeval(cv{2},basis{2},P,2)];
d2V(P>x(2,1),1)=d2V(P>x(2,1),2);
d2V(P<x(1,2),2)=d2V(P<x(1,2),1);

e=r*V-(mu*P(:,[1 1])).*dV-0.5*sigma^2*P(:,[1 1]).^2.*d2V;
e(:,2)=e(:,2)-(P-C);
e(P>x(2,1),1)=nan;
e(P<x(1,2),2)=nan;

% Plot results
vstar=[funeval(cv{1},basis{1},pstar(1));
  funeval(cv{2},basis{2},pstar(2))];
figure
plot(P,V)
hold on; plot(pstar(1),vstar(1),'*',pstar(2),vstar(2),'*'); hold off
title('Value Functions')
xlabel('P')
xlim([0 3])
ylabel('V(P)')
legend('Inactive','Active',4)

dvstar=[funeval(cv{1},basis{1},pstar(1),1);
  funeval(cv{2},basis{2},pstar(2),1)];
figure
dV(P>pstar(1),1)=NaN;
dV(P<pstar(2),2)=NaN;
plot(P,dV);
hold on; plot(pstar(1),dvstar(1),'*',pstar(2),dvstar(2),'*'); hold off
title('Marginal Value Functions')
xlabel('P')
ylabel('V''(P)')
xlim([0 3])
legend('Inactive','Active',4)

figure
plot(P,e)
title('Approximation Residual')
xlabel('P')
ylabel('Residual')
ylim([-8e-6 8e-6])
legend('Inactive','Active',4)

% Compute "closed form" solutions
optset('broyden','showiters',0)
optset('broyden','maxit',100)
[residual,beta1,beta2,A1,A2,Pl,Ph] = entexgbm([],r,mu,sigma,I,E,C);

disp('Boundary points: "Exact", approximate and error')
disp([[Ph;Pl] pstar [Ph;Pl]-pstar])

% Save Plots as EPS Files
printfigures(mfilename,3)


function out=func(flag,s,x,r,mu,sigma,C,E,I)
switch flag
  case 'f'
    out=(s-C).*(x==2);
  case 'g'
    out=mu*s;
  case 'sigma'
    out=sigma*s;
  case 'rho'
    out=r+zeros(size(s,1),1);
  case 'reward'
    out=[0 0 0;-I 0 0;-E 0 0;0 0 0];
  otherwise
    error('invalid flag')
end


% ENTEXGBM Finds entry and exit triggers for Geometric Brownian Motion
% ref.: Dixit & Pindyck, p. 218
function [residual,beta1,beta2,A1,A2,pl,ph] = ...
  entexgbm(p,rho,mu,sigma,I,E,c,pfactor,zl,zh,beta1,beta2)

if nargin==12   % ******** Computes residuals for ROOT ***********
  pl=exp(p(1,:));  % converted to logs to avoid negative values
  ph=exp(p(2,:));
  A1=p(3,:);
  A2=p(4,:);
  
  % Compute residuals from value matching & smooth pasting conditions
  residual=[
    (A1.*ph.^beta1 - A2.*ph.^beta2 - ph*pfactor + zh);
    (A1.*pl.^beta1 - A2.*pl.^beta2 - pl*pfactor + zl);
    (beta1*A1.*ph.^(beta1-1) - beta2*A2.*ph.^(beta2-1) - pfactor);
    (beta1*A1.*pl.^(beta1-1) - beta2*A2.*pl.^(beta2-1) - pfactor)];
  
elseif nargin==7 % ******** Set up problem & call BROYDEN ********
  
  if mu>=rho || sigma<=0 || c<0 || rho<=0 || I<0
    error('Improper parameter values');
  end
  
  pfactor=1/(rho-mu);
  zh=c/rho+I;
  zl=c/rho-E;
  
  s2=sigma.^2;
  beta1=0.5-mu/s2+sqrt((0.5*s2-mu).^2+2*rho*s2)/s2;
  beta2=0.5-mu/s2-sqrt((0.5*s2-mu).^2+2*rho*s2)/s2;
  
  % beta=roots([s2/2 mu-s2/2 -rho]);
  
  if isempty(p), p=[0.1;1;1;1]; end   % initialize P if starting value empty
  
  mfilename
  p=broyden(@entexgbm,log(p),rho,mu,sigma,I,E,c,pfactor,zl,zh,beta1,beta2);
  % compute residuals - should be close to zero
  residual=entexgbm(p,rho,mu,sigma,I,E,c,pfactor,zl,zh,beta1,beta2);
  
  pl=exp(p(1,:));  % converted to logs to avoid negative values
  ph=exp(p(2,:));
  A1=p(3,:);
  A2=p(4,:);
  
else
  error('Wrong number of input arguments')
end