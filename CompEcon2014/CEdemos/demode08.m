%% DEMODE08 Renewable Resource Model
%
%  Solve
%    s' = g(s)-q
%    q' = ((rho-g'(s))*(p(q)-c(s)))/p'(q)]
%  where
%    s: resource stock
%    q: harvest rate
%    x = [s;q]

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


%% STEADY-STATE

% Compute Steady State
sst = (alpha-rho)/(2*alpha);
qst = alpha*sst.*(1-sst);
xst = [sst;qst];
disp('Steady State')
disp(xst)

% Check Stability
disp('Eigenvalues')
disp(eig(fdjac(f,xst)))


%% PHASE DIAGRAM

% Layout
slim = [0.0 1.0];
qlim = [0.0 0.6];
figlayout(slim,qlim, ...
  'Renewable Resource Model Phase Diagram', ...
  'Resource Stock', ...
  'Harvest')

% Separatrix
T = 40;
xs = odespx(f,xst,T,slim,qlim);
plot(xs(:,1),xs(:,2),'k--','LineWidth',2)

% Nullclines
s = nodeunif(100,slim(1),slim(2));
snull = g(s);
plot(s,snull,[sst sst],qlim,'LineWidth',2)
legend('Separatrix','Stock Nullcline','Harvest Nullcline')
set(legend,'Box','off')
set(legend,'Location','NorthEast')
set(legend,'FontSize',12)

% Velocity Field
odefield(f,slim,qlim)


%% SOLVE ODE USING ODERK4

% Set time discretization
N  = 100;           % number of time nodes

% Find initial state on separatrix
s0 =  0.2;          % initial resource stock
q0sep = interp1(xs(:,1),xs(:,2),s0);

% Solve for different initial values and plot solution
for iv=1:3
  % Set initial values and horizon
  switch iv
    case 1          % initial value above separatrix
      T  = 2;                                % time horizon
      q0 = q0sep+0.05;                       % initial harvest
    case 2          % initial value below separatrix
      T  = 6;                                % time horizon
      q0 = q0sep-0.05;                       % initial harvest
    case 3          % initial value on separatrix
      T  = 6;                                % time horizon
      q0 = q0sep;                            % initial harvest
  end
  % Solve ODE
  x = oderk4(f,[s0;q0],T,N);
  % Plot state path
  text(x(1,1)-0.03,x(1,2),['\bf' num2str(iv)])
  bullet(x(1,1),x(1,2)+0.001,15,'r')
  odeplot(x,slim,qlim)
end
bullet(sst,qst+0.001,20)


%% SOLVE ODE USING ODECOL

% Solve ODE
n = 20;      % number of basis functions
c = zeros(n,2); c(1,:) = [s0 q0];
[x,t,res] = odecol(f,[s0;q0],T,n);

% Plot Solution
figure
plot(t,x)
title('Renewable Resource Model')
xlabel('Time')
ylabel('Quantity')
legend('Resource Stock','Harvest')
set(legend,'Box','off')
set(legend,'Location','NorthWest')
set(legend,'FontSize',14)
xlim([0 T])

% Plot Residual
figure
plot(t,res,t,0*res,'k:')
title('Renewable Resource Model')
xlabel('Time')
ylabel('Residual')
legend('Resource Stock','Harvest')
set(legend,'Box','off')
set(legend,'Location','NorthWest')
set(legend,'FontSize',14)
xlim([0 T])


%% Save Plots as EPS Files
printfigures(mfilename,3)