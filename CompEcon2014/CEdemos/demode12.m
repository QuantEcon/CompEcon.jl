%% DEMODE12 Regional Migration Model
%
%  Solve
%    L' = L*sqrt(max(0,B-k0)/alpha)
%    B' = rho*B-L.^(-theta)+wbar
%  where
%    L: labor
%    B: benefit from migration
%    x = [L;B]

% Preliminary tasks
close all
clear all


%% FORMULATION

% Model Parameters
theta = 0.5;        % inverse labor demand elasticity
k0    = 0.1;        % minimum migration cost
k1    = 0.1;        % migration cost parameter
k2    = 0.1;        % migration cost parameter
rho   = 0.05;       % discount rate
wbar  = 1;          % world wage rate

% Velocity Function
L = @(x) x(1,:);
B = @(x) x(2,:);
f = @(x) [L(x).*(sqrt(k1^2+4*k2*(B(x)-k0))-k1)/(2*k2); rho*B(x)-L(x).^(-theta)+wbar];


%% STEADY-STATE

% Compute Steady State
Bst = k0;
Lst = (wbar+rho*Bst)^(-1/theta);
xst = [Lst;Bst];
disp('Steady State')
disp(xst)

% Check Stability
disp('Eigenvalues')
disp(eig(fdjac(f,xst)))

%% PHASE DIAGRAM

% Layout
Lmin = 0.6;
Lmax = 1.4;
Bmin = 0.05;
Bmax = 0.20;
Llim = [Lmin Lmax];
Blim = [Bmin Bmax];
figlayout(Llim,Blim, ...
  'Regional Migration Model Phase Diagram', ...
  'Labor', ...
  'Benefit')

% Separatrix
T = 10;
xs = odespx(f,xst,T,Llim,Blim);
i = find(xs(:,2)>Bst);
xs = xs(i,:);
plot(xs(:,1),xs(:,2),'k--','LineWidth',2)

% Nullclines
L = nodeunif(100,Llim(1),Llim(2));
Bnull = (L.^(-theta)-wbar)/rho;
plot(L,Bnull,Llim,[Bst Bst],'LineWidth',2)
legend('Separatrix','Labor Nullcline','Benefit Nullcline')
set(legend,'Box','off')
set(legend,'Location','SouthEast')
set(legend,'FontSize',12)

% Velocity Field
odefield(f,Llim,Blim)


%% SOLVE ODE USING ODERK4

% Set time discretization
N  = 500;           % number of time nodes

% Find initial state on separatrix
L0    = 0.8;                            % initial labor
B0sep = interp1(xs(:,1),xs(:,2),L0);    % initial benefit

% Solve for different initial values and plot solution
for iv=1:3
  % Set initial values and horizon
  switch iv
    case 1          % initial value above separatrix
      T  = 1.5;                             % time horizon
      B0 = B0sep+0.02                      % initial benefit
    case 2          % initial value below separatrix
      T  = 0.35;                             % time horizon
      B0 = B0sep-0.02                      % initial benefit
    case 3          % initial value on separatrix
      T  = 2;                             % time horizon
      B0 = B0sep                           % initial benefit
  end
  % Solve ODE
  x = oderk4(f,[L0;B0],T,N);
  % Plot state path
  text(x(1,1)-0.03,x(1,2),['\bf' num2str(iv)])
  bullet(x(1,1),x(1,2),15,'r')
  odeplot(x,Llim,Blim)
end
bullet(Lst,Bst,20)


%% Save Plots as EPS Files
printfigures(mfilename,1)