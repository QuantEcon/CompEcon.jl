%% DEMODE04 - Non-IVP Non-Homogeneous Linear ODE Example
%
%  Solve x1' = -1*x1 - 0.5*x2 + 2
%        x2' =        -0.5*x2 + 1
%  s.t   x1(0)=1, x2(1)=1
%        t in [0,10]

% Preliminary tasks
demosetup(mfilename)


%% FORMULATION

% Model Parameters
A = [-1 -0.5; 0 -0.5];      % velocity function parameters
b = [2; 1];                 % velocity function parameters

% Velocity Function
f = @(x) [A(1,1)*x(1,:)+A(1,2)*x(2,:)+b(1); ...
          A(2,1)*x(1,:)+A(2,2)*x(2,:)+b(2)];
        
% Closed-Form Solution 
X = @(t) [1-exp(1/2-t)+exp((1-t)/2) 2-exp((1-t)/2)];
                

%% STEADY-STATE

% Compute Steady State
xst = -A\b;
disp('Steady State')
disp(xst)
disp('Eigenvalues')
disp(eig(A))


%% SOLVE ODE ANALYTICALLY

% Plot Solution
figure
t = nodeunif(200,0,T);
plot(t,X(t))
title('Solution')
xlabel('Time')
legend('x_1','x_2',2)
set(legend,'Box','off')
set(legend,'Location','NorthWest')
set(legend,'FontSize',14)


%% SOLVE ODE USING ODECOL

% Solve ODE
n = 15;                       % number of basis functions
T = 10;                       % time horizon
bx = [1;2];                   % boundary variables
bt = [0;1];                   % boundary times
bv = [1;1];                   % boundary values
[x,t,res] = odecol(f,bv,T,n,bt,bx);

% Plot Appproximation Errors
figure
plot(t,x-X(t),t,0*t,'k:')
title('Approximation Errors')
xlabel('Time')
ylabel('Error')
legend('x_1''','x_2''')
set(legend,'Box','off')
set(legend,'Location','NorthWest')
set(legend,'FontSize',14)

% Plot Residuals
figure
plot(t,res,t,0*res,'k:')
title('Approximation Residuals')
xlabel('Time')
ylabel('Residual')
legend('x_1''','x_2''')
set(legend,'Box','off')
set(legend,'Location','NorthWest')
set(legend,'FontSize',14)


%% PHASE DIAGRAM

% Layout
x1lim = [0 2];
x2lim = [0 4];
figlayout(x1lim,x2lim,'Phase Diagram','x_1','x_2')

% Nullclines
x1 = nodeunif(100,x1lim(1),x1lim(2));
x1null = -(A(1,1)*x1+b(1))/A(1,2); 
x2null = -(A(2,1)*x1+b(2))/A(2,2); 
plot(x1,x1null,x1,x2null,'LineWidth',2)
legend('x_1 Nullcline','x_2 Nullcline')
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
printfigures(mfilename,4)