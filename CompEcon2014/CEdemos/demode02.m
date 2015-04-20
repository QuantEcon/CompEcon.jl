%% DEMODE02 Generic IVP Nonlinear ODE Example
%
%  Solve
%    x1' = x1^2 - 2x2 - a
%    x2' = b - x1 - x2

% Preliminary tasks
demosetup(mfilename)


%% FORMULATION

% Model Parameters
a =  5;
b = -1;

% Velocity Function
f = @(x) [x(1,:).^2-2*x(2,:)-a; b-x(1,:)-x(2,:)];


%% STEADY STATES

% Steady State A
A = zeros(2,1);
A(1) = (-2-sqrt(4+4*(2*b+a)))/2;
A(2) = b-A(1);
disp('Steady State A')
disp(A)
J = [2*A(1) -2;-1 -1];
disp('Eigenvalues')
disp(eig(J))
0.5*(-7+sqrt(33))
0.5*(-7-sqrt(33))

% Steady State B
B = zeros(2,1);
B(1) = (-2+sqrt(4+4*(2*b+a)))/2;
B(2) = b-B(1);
disp('Steady State B')
disp(B)
J = [2*B(1) -2;-1 -1];
disp('Eigenvalues')
disp(eig(J))


%% PHASE DIAGRAM

% Layout
x1lim = [A(1)-2 B(1)+3];
x2lim = [B(2)-4 A(2)+2];
figlayout(x1lim,x2lim,'Nonlinear ODE Phase Diagram','x_1','x_2')

% Separatrix
T = 8;
xs = odespx(f,B,T,x1lim,x2lim);
plot(xs(:,1),xs(:,2),'k--','LineWidth',2)

% Nullclines
x1 = nodeunif(100,x1lim(1),x1lim(2));
x2 = [(x1.^2-5)/2 -1-x1];
plot(x1,x2,'LineWidth',2)

% Velocity Field
odefield(f,x1lim,x2lim)


%% SOLVE ODE USING ODERK4

% Set time discretization
N  = 100;           % number of time nodes

% Find initial state on separatrix
x1sep = -1;
x2sep = interp1(xs(:,1),xs(:,2),x1sep);

% Solve for different initial values and plot solution
for iv=1:3
  % Set initial values and horizon
  switch iv
    case 1          % initial value above separatrix
      T  = 8;                         % time horizon
      x0 = [x1sep;x2sep+0.7];         % initial values
    case 2          % initial value below separatrix
      T  = 4;                         % time horizon
      x0 = [x1sep;x2sep-0.7];         % initial values
    case 3          % initial value on separatrix
      T  = 3;                         % time horizon
      x0 = [x1sep;x2sep];             % initial values
  end
  % Solve ODE
  x = oderk4(f,x0,T,N);
  % Plot state path
  text(x(1,1)-0.3,x(1,2),['\bf' num2str(iv)])
  bullet(x(1,1),x(1,2)+0.05,20,'r')
  odeplot(x,x1lim,x2lim)
end

% Plot Steady States
bullet(A(1),A(2),20)
text(A(1)+0.2,A(2)+0.1,'A','FontSize',14)
bullet(B(1),B(2),20)
text(B(1)+0.4,B(2)+0.1,'B','FontSize',14)


%% Save Plots as EPS Files
printfigures(mfilename,1)