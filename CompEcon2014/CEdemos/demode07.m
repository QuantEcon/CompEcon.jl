%% DEMODE07 Commercial Fisheries Model (V.L. Smith) 
%  Intertemporal price equilibrium in market for renewable resource.
%
%  Solve
%    s' = (1-s)*s - s*k/(alpha+beta*s*k) 
%    k' = delta*(alpha*s/(2(alpha+beta*s*k)^2}-phi)
%  where
%    s: Stock of fish
%    k: Industry size
%    x = [s;k]

% Preliminary tasks
demosetup(mfilename)


%% FORMULATION

% Model Parameters
alpha = 0.5;    % marginal cost parameter
beta  = 2.75;   % slope of demand curve
phi   = 0.05;   % fixed cost
delta = 10;     % industry entry rate

% Velocity Function
s = @(x) x(1,:);
k = @(x) x(2,:);
y = @(x) s(x)./(alpha+beta*s(x).*k(x));
p = @(x) alpha./(alpha+beta*s(x).*k(x));
f = @(x) [(1-s(x)).*s(x)-k(x).*y(x); delta*(0.5*p(x).*y(x)-phi)];


%% STEADY-STATE

% Steady State A
A = [0.1;0.6]; 
A = broyden(f,A);
disp('Steady State A')
disp(A)
J = fdjac(f,A);
disp('Eigenvalues')
disp(eig(J))

% Steady State B
B = [0.3;0.9]; 
B = broyden(f,B);
disp('Steady State B')
disp(B)
J = fdjac(f,B);
disp('Eigenvalues')
disp(eig(J))

% Steady State C
C = [0.5;0.8]; 
C = broyden(f,C);
disp('Steady State C')
disp(C)
J = fdjac(f,C);
disp('Eigenvalues')
disp(eig(J))


%% PHASE DIAGRAM

% Layout
slim = [0 1];
klim = [0 1];
figlayout(slim,klim, ...
  'Commercial Fisheries Model Phase Diagram', ...
  'Fish Population', ...
  'Industry Size')

% Separatrix
T = 12;
xs = odespx(f,B,T,slim,klim);
plot(xs(:,1),xs(:,2),'k--','LineWidth',2)

% Nullclines
s = nodeunif(100,slim(1),slim(2));
snull = alpha*(1-s)./(1-beta*s.*(1-s)); 
knull = (sqrt(alpha*s/(2*phi))-alpha)./(beta*s);
plot(s,snull,s,knull,'LineWidth',2)
legend('Separatrix','Fish Nullcline','Industry Nullcline')
set(legend,'Box','on')
set(legend,'Location','South')
set(legend,'FontSize',12)

% Velocity Field
odefield(f,slim,klim)


%% SOLVE ODE USING ODERK4

% Set time discretization
N  = 200;           % number of time nodes

% % Solve for different initial values and plot solution
% for iv=1:3
%   % Set initial values and horizon
%   switch iv
%     case 1          % initial value left of separatrix
%       T  = 40;                              % time horizon
%       s0 = 0.2;                             % initial fish population
%       k0 = 0.0;                             % initial industry size
%     case 2          % initial value on separatrix
%       T  = 5;                               % time horizon
%       s0 = 0.3;                             % initial fish population
%       k0 = interp1(xs(:,1),xs(:,2),s0);     % initial industry size
%     case 3          % initial value right of separatrix
%       T  = 18;                              % time horizon
%       s0 = 0.8;                             % initial fish population
%       k0 = 0.0;                             % initial industry size
%   end
%   % Solve ODE
%   x = oderk4(f,[s0;k0],T,N);
%   % Plot state path
%   bullet(x(1,1),x(1,2)+0.001,15,'r')
%   odeplot(x,slim,klim)
% end

% Solve for different initial values and plot solution
N  = 400;
T  =  30;
Ns =  11;
Nk =   2;
sinit = nodeunif(Ns,slim(1),slim(2)); sinit(1)=[];
kinit = nodeunif(Nk,klim(1),klim(2));
for is=1:Ns-1
  for ik=1:Nk
    x = oderk4(f,[sinit(is);kinit(ik)],T,N);
    bullet(x(1,1),x(1,2)+0.001,15,'r')
    plot(x(:,1),x(:,2),'m','LineWidth',2)
  end
end

% Plot steady states
bullet(A(1),A(2)+0.01,20)
text(A(1)+0.03,A(2)-0.01,'A','FontSize',14)
bullet(B(1),B(2)+0.01,20)
text(B(1)+0.01,B(2)-0.07,'B','FontSize',14)
bullet(C(1),C(2)+0.01,20)
text(C(1)-0.04,C(2)-0.06,'C','FontSize',14)


%% SOLVE ODE USING ODECOL

% Solve ODE
T  = 40;                              % time horizon
s0 = 0.2;                             % initial fish population
k0 = 0.0;                             % initial industry size
n = 100;                               % number of basis functions
% c = zeros(n,2); c(1,:) = B;
% [x,t,res] = odecol(f,[s0;k0],T,n,[],[],[],c);
c = zeros(n,2); c(1,:) = [s0 .1];
[x,t,res] = odecol(f,[s0;k0],T,n,[],[],[],c);

% Plot Solution
figure
plot(t,x)
title('Commercial Fisheries Model')
xlabel('Time')
ylabel('Population/Size')
legend('Fish Population','Industry Size')
set(legend,'Box','off')
set(legend,'Location','NorthEast')
set(legend,'FontSize',14)

% Plot Residual
figure
plot(t,res,t,0*res,'k:')
title('Commercial Fisheries Model')
xlabel('Time')
ylabel('Residual')
legend('Fish Population','Industry Size')
set(legend,'Box','off')
set(legend,'Location','NorthWest')
set(legend,'FontSize',14)


%% Save Plots as EPS Files
printfigures(mfilename,3)