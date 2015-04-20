%% DEMODE02HW Generic Nonlinear ODE Example Variant for Homework Assignment
%
%  Solve
%    x1' = x1^2 + 2x2 - a
%    x2' = b - x1 - x2
%  s.t.
%    x(0) given

% Preliminary tasks
demosetup(mfilename)


%% FORMULATION

% Model Parameters
a = 5;
b = 1;

% Velocity Function
f = @(x) [x(1,:).^2+2*x(2,:)-a; b-x(1,:)-x(2,:)];


%% STEADY STATES

% Steady State A
A = zeros(2,1);
A(1) = (2+sqrt(4-4*(2*b-a)))/2;
A(2) = b-A(1);
disp('Steady State A')
disp(A)
J = [2*A(1) 1;-1 -1];
disp('Eigenvalues')
disp(eig(J))

% Steady State B
B = zeros(2,1);
B(1) = (2-sqrt(4-4*(2*b-a)))/2;
B(2) = b-B(1);
disp('Steady State B')
disp(B)
J = [2*B(1) 1;-1 -1];
disp('Eigenvalues')
disp(eig(J))


%% PHASE DIAGRAM

% Layout
x1lim = [B(1)-2 A(1)+2];
x2lim = [A(2)-2 B(2)+2];
figlayout(x1lim,x2lim,'Nonlinear ODE Phase Diagram','x_1','x_2')

% Separatrix
T = 20;
xs = odespx(f,A,T,x1lim,x2lim);
plot(xs(:,1),xs(:,2),'k--','LineWidth',2)

% Nullclines
x1node = nodeunif(100,x1lim(1),x1lim(2));
x1null = (a-x1node.^2)/2;
x2null = b - x1node;
plot(x1node,x2null,'g','LineWidth',2)
plot(x1node,x1null,'b','LineWidth',2)

% Velocity Field
odefield(f,x1lim,x2lim)


%% SOLVE ODE USING ODERK4

% Set time discretization
N  = 200;           % number of time nodes

% Solve for different initial values and plot solution
for iv=1:3
  % Set initial values and horizon
  switch iv
    case 1          % initial value left of separatrix
      T  = 5;                         % time horizon
      x0 = [-2; -2];                  % initial values
    case 2          % initial value left of separatrix
      T  = 3;                         % time horizon
      x0 = [-2; 3];                    % initial values
    case 3          % initial value right of separatrix
      T  = 0.5;                       % time horizon
      x0 = [3; 2];                    % initial values
  end
  % Solve ODE
  x = oderk4(f,x0,T,N);
  % Plot state path
  text(x(1,1)-0.2,x(1,2)+0.2,['\bf' num2str(iv)])
  bullet(x(1,1),x(1,2)+0.05,20,'r')
  odeplot(x,x1lim,x2lim)
end

% Plot Steady States
bullet(A(1),A(2),20)
text(A(1)+0.2,A(2)+0.1,'A','FontSize',14)
bullet(B(1),B(2),20)
text(B(1)+0.4,B(2)+0.1,'B','FontSize',14)


%% Save Plots as EPS Files
% printfigures(mfilename,4)