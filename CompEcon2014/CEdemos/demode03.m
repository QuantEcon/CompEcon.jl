%% DEMODE03 - Initial Value Non-Homogeneous Linear ODE Example
%
%  Solve x1' =  1*x1 + 12*x2 - 60
%        x2' = -1*x1 -  6*x2 + 36
%  s.t.  x1(0)=5, x2(0)=2
%        t in [0,3]

% Preliminary tasks
demosetup(mfilename)


%% FORMULATION

% Model Parameters
A  = [1 12; -1 -6];     % velocity function parameters
b  = [-60; 36];         % velocity function parameters

% Velocity Function
f = @(x) [A(1,1)*x(1,:)+A(1,2)*x(2,:)+b(1); ...
          A(2,1)*x(1,:)+A(2,2)*x(2,:)+b(2)];
        
% Closed-Form Solution 
X = @(t) [12 - 48*exp(-2*t) + 42*exp(-3*t) ...
           4 + 12*exp(-2*t) - 14*exp(-3*t)];
   
         
%% STEADY-STATE

% Compute Steady State
xst = -A\b;
disp('Steady State')
disp(xst)
disp('Eigenvalues')
disp(eig(A))

                  
%% CLOSED FORM SOLUTION

% Plot Solution
figure
N = 200;                % number of time nodes
T = 3;                  % time horizon
t = nodeunif(N,0,T);
plot(t,X(t))
title('Solution')
xlabel('Time')
legend('x_1','x_2',2)
set(legend,'Box','off')
set(legend,'Location','NorthWest')
set(legend,'FontSize',14)


%% SOLVE ODE USING ODERK4

% Solve ODE
tic
N  = 1000;              % number of time nodes
T  = 3;                 % time horizon
x0 = [6;2];             % initial values
[x,t] = oderk4(f,x0,T,N);
toc

% Plot Approximation Errors
figure
plot(t,x-X(t),t,0,'k:')
title('Approximation Errors - Runge Kutta')
xlabel('Time')
ylabel('Error')
legend('x_1','x_2',0)
set(legend,'Box','off')
set(legend,'Location','NorthEast')
set(legend,'FontSize',14)


%% SOLVE ODE USING ODECOL

% Solve ODE
optset('odecol','showiters',1)
tic
n  = 20;                % number of basis functions
T  = 3;                 % time horizon
x0 = [6;2];             % initial values
[x,t,res] = odecol(f,x0,T,n);
toc

% Plot Appproximation Errors
figure
plot(t,x-X(t),t,0*t,'k:')
title('Approximation Errors - Collocation')
xlabel('Time')
ylabel('Error')
legend('x_1','x_2')
set(legend,'Box','off')
set(legend,'Location','NorthEast')
set(legend,'FontSize',14)

% Plot Residuals
figure
plot(t,res,t,0*res,'k:')
title('Approximation Residuals - Collocation')
xlabel('Time')
ylabel('Residual')
legend('x_1','x_2')
set(legend,'Box','off')
set(legend,'Location','NorthEast')
set(legend,'FontSize',14)


%% PHASE DIAGRAM

% Layout
x1lim = [0 15];
x2lim = [0  8];
figlayout(x1lim,x2lim,'Phase Diagram','x_1','x_2')

% Separatrix
T = 12;
xs = odespx(f,xst,T,x1lim,x2lim);
plot(xs(:,1),xs(:,2),'k--','LineWidth',2)

% Nullclines
x1 = nodeunif(100,x1lim(1),x1lim(2));
x1null = -(A(1,1)*x1+b(1))/A(1,2); 
x2null = -(A(2,1)*x1+b(2))/A(2,2); 
plot(x1,x1null,x1,x2null,'LineWidth',2)
legend('Separatrix','x_1 Nullcline','x_2 Nullcline')
set(legend,'Box','off')
set(legend,'Location','NorthEast')
set(legend,'FontSize',14)

% State Path
bullet(x(1,1),x(1,2)+0.02,20,'r')
odeplot(x,x1lim,x2lim)
bullet(xst(1),xst(2)+0.02,20)

% Velocity Field
odefield(f,x1lim,x2lim)


%% Save Plots as EPS Files
printfigures(mfilename,5)