%% DEMODE09 Economic Growth Model
%
%  Solve
%    k' = k^alpha-delta*k-q
%    q' = (alpha*k^{alpha-1}-delta-rho)*q/theta
%  where
%    k: capital stock
%    q: consumption
%    x = [k;q]

% Preliminary tasks
demosetup(mfilename)


%% FORMULATION

% Model Parameters
alpha = 0.5;               % capital share
delta = 0.05;              % depreciation rate
theta = 2.0;               % relative risk aversion
rho   = 0.05;              % discount rate

% Velocity Function
f = @(x) [x(1,:).^alpha-delta*x(1,:)-x(2,:); ...
         (alpha*x(1,:).^(alpha-1)-delta-rho).*x(2,:)/theta];


%% STEADY-STATE

% Compute Steady State
kst = ((delta+rho)/alpha)^(1/(alpha-1));
qst = kst^alpha-delta*kst;
pst = qst^-theta;
vst = ((qst.^(1-theta))/(1-theta))/rho;
xst = [kst;qst];
disp('Steady State')
disp(xst)

% Check Stability
disp('Eigenvalues')
disp(eig(fdjac(f,xst)))


%% PHASE DIAGRAM

% Layout
kmin = 10;
kmax = 30;
qmin = 2;
qmax = 4;
klim = [kmin kmax];
qlim = [qmin qmax];
figlayout(klim,qlim, ...
  'Economic Growth Model Phase Diagram', ...
  'Capital Stock', ...
  'Consumption')

% Separatrix
T = 300;
xs = odespx(f,xst,T,klim,qlim);
plot(xs(:,1),xs(:,2),'k--','LineWidth',2)

% Nullclines
k = nodeunif(100,klim(1),klim(2));
qnull = k.^alpha-delta*k;
plot(k,qnull,[kst kst],qlim,'LineWidth',2)
legend('Separatrix','Capital Nullcline','Consumption Nullcline')
set(legend,'Box','off')
set(legend,'Location','NorthWest')
set(legend,'FontSize',12)

% Velocity Field
odefield(f,klim,qlim)


%% SOLVE ODE USING ODERK4

% Set time discretization
N  = 100;           % number of time nodes

% Find initial state on separatrix
k0    = 15;         % initial capital stock
q0sep = interp1(xs(:,1),xs(:,2),k0);

% Solve for different initial values and plot solution
for iv=1:3
  % Set initial values and horizon
  switch iv
    case 1          % initial value above separatrix
      T  = 12;                              % time horizon
      q0 = q0sep+0.2;                           % initial consumption
    case 2          % initial value below separatrix
      T  = 20;                              % time horizon
      q0 = q0sep-0.2;                           % initial consumption
    case 3          % initial value on separatrix
      T  = 90;                               % time horizon
      q0 = q0sep;                           % initial consumption
  end
  % Solve ODE
  x = oderk4(f,[k0;q0],T,N);
  % Plot state path
  text(x(1,1)-0.7,x(1,2),['\bf' num2str(iv)])
  bullet(x(1,1),x(1,2),15,'r')
  odeplot(x,klim,qlim)
end
bullet(kst,qst+0.001,20)


%% SOLVE ODE USING ODECOL

% Solve ODE
n = 20;                            % number of basis functions
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
xlim([0 T])

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
xlim([0 T])


%% Save Plots as EPS Files
printfigures(mfilename,3)