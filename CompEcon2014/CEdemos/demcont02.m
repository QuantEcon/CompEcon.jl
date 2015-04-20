function demcont02

%% DEMCONT02 Stochastic Renewable Resource Model
%
% Old demode08alt and demsc02

% Preliminary tasks
demosetup(mfilename)


%% FORMULATION

% Model Parameters
alpha = 0.5;        % biological growth function scale factor
beta  = 0.5;        % biological growth function elasticity factor
kappa = 5;          % unit harvest cost scale factor
gamma = 1.5;        % unit harvest cost elasticity
eta   = 1.5;        % inverse elasticity of demand
rho   = 0.05;       % discount rate
sigma = 0.1;        % diffusion volatiity

% Ancillary Functions
s    = @(x) x(1,:);                             % map x to s
q    = @(x) x(2,:);                             % map x to q
p    = @(q) q.^(-eta);                          % inverse demand func
pder = @(q) -eta*q.^(-eta-1);                   % inverse demand deriv
g    = @(s) (alpha/beta)*s.*(1-s.^beta);        % biological growth func
gder = @(s) (alpha/beta)*(1-(1+beta)*s.^beta);  % biological growth deriv
k    = @(s) kappa*s.^(-gamma);                  % unit harvest func
kder = @(s) -kappa*gamma*s.^(-gamma-1);         % unit harvest deriv

% Velocity Function
f = @(x) [g(s(x))-q(x); ...
   ((rho-gder(s(x))).*(p(q(x))-k(s(x)))+kder(s(x)).*g(s(x)))./pder(q(x))];
 

%% DETERMINISTIC STEADY-STATE

% Compute Steady State
sst = 0.6;
qst = g(sst);
xst = [sst;qst];
xst = broyden(f,xst);
disp('Steady State')
disp(xst)
sst = xst(1);
qst = xst(2);

% Check Stability
disp('Eigenvalues')
disp(eig(fdjac(f,xst)))


%% SOLVE USING SCSOLVE

% Approximation Structure
n = 1000;                            	% number of basis functions
smin = 0.1;                             % minimum resource stock
smax = 1.0;                             % maximum resource stock
basis = fundefn('spli',n,smin,smax);    % basis functions

% Model Structure
model.func   = @func;                   % model function file
model.rho    = rho;                     % discount rate
model.params = {g,k,eta,sigma};         % model parameters

% Solve Bellman Equation
[c,s,v,q,resid] = scsolve(model,basis);

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


%% ERGODIC DISTRIBUTION

figure
cp = itodensity(model,basis,c);
plot(s,max(0,funeval(cp,basis,s)))
title('Ergodic Distribution')
xlabel('Resource Stock')
ylabel('Probability')
xlim([smin smax])


%% PLOT ERRORS FOR CASES WITH ANALYTIC SOLUTION

% Plot Percent Approximation Error for Cases with Closed-Form Solution
if gamma==(1+beta)&&eta==(1+beta)
  optset('broyden','tol',1e-14)
  optset('broyden','showiters',1)
  theta = ((alpha+rho)*beta/(1+beta)-0.5*beta*sigma^2)^((1+beta)/beta);
  fphi = @(phi) theta*phi^((1+beta)/beta)-beta*phi-kappa;
  optset('broyden','tol',1e-12)
  phi = broyden(fphi,12);
  va = -phi*(s.^-beta+alpha/rho);
  qa = ((kappa+beta*phi)^(-1/(1+beta)))*s;
  figure
  hold on
  plot(s,v./va-1,s,q./qa-1)
  legend('Value Function','Optimal Policy','Location','Best')
  legend boxoff
  plot(s,0*s,'k--','Linewidth',1)
  title('Percent Approximation Errors')
  xlabel('Resource Stock')
  ylabel('Error')
end


%% Save Plots as EPS Files
printfigures(mfilename,5)


%% Function File for scsolve
function out = func(flag,s,q,Vs,g,k,eta,sigma)
switch flag
  case 'x'
    out = (Vs+k(s)).^(-1/eta);
  case 'f'
    out = (1/(1-eta))*q.^(1-eta) - k(s).*q;
  case 'g'
    out = g(s) - q;
  case 'sigma'
    out = sigma*s;
end