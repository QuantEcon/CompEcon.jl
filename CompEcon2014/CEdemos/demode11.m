%% DEMODE11 Tobin's Q Model
%
%  Solve
%    k' = k.*(max(q-1,0)/b-delta)
%    q' = (delta+rho)*q-1./k-(q-1)^2/(2*b)]
%  where
%    k: capital stock
%    q: price of capital
%    x = [k;q]

% Preliminary tasks
demosetup(mfilename)


%% FORMULATION

% Model Parameters
delta = 0.05;              % depreciation rate
b     = 5;                 % marginal adjustment cost
rho   = 0.05;              % interest rate

% Velocity Function
k = @(x) x(1,:);
q = @(x) x(2,:);
f = @(x) [k(x).*(max(q(x)-1,0)/b-delta); ...
  (delta+rho)*q(x)-1./k(x)-(q(x)-1).^2/(2*b)]


%% STEADY-STATE

% Compute Steady State
qst = 1+b*delta;
kst = 1/((delta+rho)*qst-(b*delta)^2/(2*b));
xst = [kst;qst];
disp('Steady State')
disp(xst)

% Check Stability
disp('Eigenvalues')
disp(eig(fdjac(f,xst)))


%% PHASE DIAGRAM

% Layout
kmin = 7;
kmax = 10;
qmin = 1.0;
qmax = 1.5;
klim = [kmin kmax];
qlim = [qmin qmax];
figlayout(klim,qlim, ...
  'Tobin''s Q Model Phase Diagram', ...
  'Capital Stock', ...
  'Price of Capital')

% Separatrix
T = 300;
xs = odespx(f,xst,T,klim,qlim);
plot(xs(:,1),xs(:,2),'k--','LineWidth',2)

% Nullclines
k = nodeunif(100,klim(1),klim(2));
A = -1/(2*b);
B = (delta+rho)+1/b;
C = 1./k+1/(2*b);
i = find(isfinite(B^2+4*A*C)==1);
qnull = (-B+sqrt(B^2+4*A*C))/(2*A);
plot(k(i),qnull(i),[kst kst],qlim,'LineWidth',2)
legend('Separatrix','Capital Nullcline','Price Nullcline')
set(legend,'Box','off')
set(legend,'Location','NorthWest')
set(legend,'FontSize',12)

bullet(kst,qst,20)

% Velocity Field
odefield(f,klim,qlim)


%% SOLVE ODE USING ODERK4

% Set time discretization
N  = 100;           % number of time nodes

% Find initial state on separatrix
k0    = 9;          % initial capital
q0sep = interp1(xs(:,1),xs(:,2),k0);

% Solve for different initial values and plot solution
for iv=1:3
  % Set initial values and horizon
  switch iv
    case 1          % initial value above separatrix
      T  = 10;                              % time horizon
      q0 = q0sep+0.07;                      % initial price
    case 2          % initial value below separatrix
      T  = 10;                              % time horizon
      q0 = q0sep-0.07;                      % initial price
    case 3          % initial value on separatrix
      T  = 20;                              % time horizon
      q0 = q0sep;                           % initial price
  end
  % Solve ODE
  x = oderk4(f,[k0;q0],T,N);
  % Plot state path
  text(x(1,1)-0.01,x(1,2)+0.02,['\bf' num2str(iv)])
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
legend('Capital','Price')
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
legend('Capital','Price')
set(legend,'Box','off')
set(legend,'Location','NorthWest')
set(legend,'FontSize',14)
xlim([0 T])


% %% Save Plots as EPS Files
% printfigures(mfilename,3)