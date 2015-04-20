%% DEMODE09HW Economic Growth Model - Homework variant that uses shadow price rather than consumption as second ODE variable.
%
%  Solve
%    k' = k^alpha-delta*k-p^-theta
%    p' = (rho+delta-alpha*k^{alpha-1})*q
%  where
%    k: capital stock
%    p: market price of consumption good
%    x = [k;p]

% Preliminary tasks
demosetup(mfilename)


%% FORMULATION

% Model Parameters
alpha = 0.4;               % capital share
delta = 0.05;              % depreciation rate
theta = 0.5;               % inverse relative risk aversion
rho   = 0.05;              % discount rate

% Velocity Function
k = @(x) x(1,:);
p = @(x) x(2,:);
f = @(x) [k(x).^alpha-delta*k(x)-p(x).^-theta; ...
         (rho+delta-alpha*k(x).^(alpha-1)).*p(x)];


%% STEADY-STATE

% Compute Steady State
kst = ((delta+rho)/alpha)^(1/(alpha-1));
pst = (kst^alpha-delta*kst)^(-1/theta);
qst = pst^(-theta);
vst = ((qst.^(1-theta))/(1-theta))/rho;
xst = [kst;pst];
f(xst)
disp('Steady State')
disp(xst)

% Check Stability
J = fdjac(f,xst);
disp('Eigenvalues')
disp(eig(J))

%% PHASE DIAGRAM

% Layout
kmin = round(kst)-3;
kmax = round(kst)+3;
pmin = round(pst*10)/10-0.3;
pmax = round(pst*10)/10+0.3;
klim = [kmin kmax];
plim = [pmin pmax];
figlayout(klim,plim, ...
  'Economic Growth Model Phase Diagram', ...
  'Capital Stock', ...
  'Market Price')

% Separatrix
T = 200;
xs = odespx(f,xst,T,klim,plim);
plot(xs(:,1),xs(:,2),'k--','LineWidth',2)

% Nullclines
k = nodeunif(100,klim(1),klim(2));
pnull = (k.^alpha-delta*k).^(-1/theta);
plot(k,pnull,[kst kst],plim,'LineWidth',2)
legend('Separatrix','Capital Nullcline','Consumption Nullcline')
set(legend,'Box','off')
set(legend,'Location','NorthWest')
set(legend,'FontSize',14)

% Velocity Field
odefield(f,klim,plim)


%% SOLVE ODE USING ODERK4

% Set time discretization
N  = 100;           % number of time nodes

% Find initial state on separatrix
k0    =  8;         % initial capital stock
q0sep = interp1(xs(:,1),xs(:,2),k0);

% Solve for different initial values and plot solution
for iv=1:3
  % Set initial values and horizon
  switch iv
    case 1          % initial value above separatrix
      T  = 30;                              % time horizon
      q0 = q0sep+0.03;                      % initial consumption
    case 2          % initial value below separatrix
      T  = 30;                              % time horizon
      q0 = q0sep-0.03;                      % initial consumption
    case 3          % initial value on separatrix
      T  = 60;                              % time horizon
      q0 = q0sep;                           % initial consumption
  end
  % Solve ODE
  x = oderk4(f,[k0;q0],T,N);
  % Plot state path
  text(x(1,1)-0.2,x(1,2)-0.01,['\bf' num2str(iv)])
  bullet(x(1,1),x(1,2),15,'r')
  odeplot(x,klim,plim)
end
bullet(kst,qst+0.001,20)


%% SOLVE ODE USING ODECOL

% Solve ODE
T  = 60;                            % time horizon
n  = 25;                            % number of basis functions
[x,t,res] = odecol(f,[k0;q0],T,n);

% Plot Solution
figure
plot(t,x)
title('Economic Growth Model')
xlabel('Time')
ylabel('Quantity')
legend('Capital Stock','Consumption')
set(legend,'Box','off')
set(legend,'Location','NorthWest')
set(legend,'FontSize',14)

% Plot Residual
figure
plot(t,res,t,0*res,'k:')
title('Economic Growth Model')
xlabel('Time')
ylabel('Residual')
legend('Capital Stock','Consumption')
set(legend,'Box','off')
set(legend,'Location','NorthWest')
set(legend,'FontSize',14)


%% Save Plots as EPS Files
printfigures(mfilename,3)