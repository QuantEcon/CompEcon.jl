function demcont08

%% DEMCONT08 Deterministic Renewable Resource Model

% Preliminary tasks
demosetup(mfilename)


%% FORMULATION

% Model Parameters
alpha = 2;          % biological growth function scale factor
kappa = 1;          % unit harvest cost scale factor
eta   = 2;          % inverse elasticity of demand
rho   = 0.05;       % discount rate

% Ancillary Functions
s    = @(x) x(1,:);                 % map x to s
q    = @(x) x(2,:);                 % map x to q
p    = @(x) q(x).^(-eta);           % inverse demand func
pder = @(x) -eta*q(x).^(-eta-1);    % inverse demand deriv
g    = @(s) alpha*s.*(1-s);         % biological growth func
gder = @(s) alpha*(1-2*s);          % biological growth deriv

% Velocity Function
f    = @(x) [g(s(x))-q(x); ...
   ((rho-gder(s(x))).*(p(x)-kappa))./pder(x)];
 

%% DETERMINISTIC STEADY-STATE

% Compute Steady State
sst = (alpha-rho)/(2*alpha);
qst = alpha*sst.*(1-sst);
xst = [sst;qst];
disp('Steady State')
disp(xst)

% Check Stability
disp('Eigenvalues')
disp(eig(fdjac(f,xst)))


%% SOLVE USING SCSOLVE

% Approximation Structure
n = 30;                            	% number of basis functions
smin = 0.1;                             % minimum resource stock
smax = 1.0;                             % maximum resource stock
basis = fundefn('cheb',n,smin,smax);     % basis functions
s = funnode(basis);                     % collocation nodes

% Model Structure
model.func   = @func;                   % model function file
model.rho    = rho;                     % discount rate
model.params = {g,kappa,eta};           % model parameters

% Solve Bellman Equation
q = (qst/sst)*s;
v = ((1/(1-eta))*q.^(1-eta)-kappa*q)/rho;
[c,s,v,q,resid] = scsolve(model,basis,v);

% Plot Value Function
figure
hold on
vst = funeval(c,basis,sst);
plot(s,v)
plot(sst,vst,'k*','LineWidth',8)
xlim([smin smax])
title('Value Function')
xlabel('Resource Stock')
ylabel('Value')

% Plot Shadow Price Function
figure
hold on
p = funeval(c,basis,s,1);
pst = funeval(c,basis,sst,1);
plot(s,p)
plot(sst,pst,'k*','LineWidth',8)
xlim([smin smax])
title('Shadow Price Function')
xlabel('Resource Stock')
xlabel('Shadow Price')

% Plot Optimal Policy
figure
hold on
plot(s,q)
plot(sst,qst,'k*','LineWidth',8)
xlim([smin smax])
title('Optimal Policy')
xlabel('Resource Stock')
ylabel('Harvest')

% Plot Residual
figure
plot(s,resid,s,0*s,'k:')
xlim([smin smax])
title('Bellman Equation Residual')
xlabel('Resource Stock')
ylabel('Residual')


%% Save Plots as EPS Files
printfigures(mfilename,4)


%% Function File for scsolve
function out = func(flag,s,q,Vs,g,kappa,eta)
switch flag
  case 'x'
    out = (Vs+kappa).^(-1/eta);
  case 'f'
    out = (1/(1-eta))*q.^(1-eta) - kappa*q;
  case 'g'
    out = g(s) - q;
  case 'sigma'
    out = [];
end