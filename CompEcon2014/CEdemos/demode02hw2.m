%% DEMODE02HW Generic Nonlinear ODE Example Variant for Homework Assignment
%
%  Solve
%    x1' = (x1-a1)*(x2-a2)
%    x2' =  b0 - b1*x1 - b2*x2

% Preliminary tasks
demosetup(mfilename)


%% FORMULATION

% Model Parameters
a1 = 1;
a2 = 2;
b0 = 14;
b1 = 2;
b2 = 4;  

% Velocity Function
f = @(x) [(x(1,:)-a1).*(x(2,:)-a2); b0-b1*x(1,:)-b2*x(2,:)];


%% STEADY STATES

% Steady State A
A = [a1;b0/b2-(b1/b2)*a1];
disp('Steady State A')
disp(A)
J = [A(2)-a2 A(1)-a1;-b1 -b2];
disp('Eigenvalues')
disp(eig(J))

% Steady State B
B = [b0/b1-(b2/b1)*a2;a2];
disp('Steady State B')
J = [B(2)-a2 B(1)-a1;-b1 -b2];
disp('Eigenvalues')
disp(eig(J))



%% PHASE DIAGRAM

% Layout
x1lim = [-2 4];
x2lim = [ 0 5];
figlayout(x1lim,x2lim,'Nonlinear ODE Phase Diagram','x_1','x_2')

% Separatrix
T = 4;
xs = odespx(f,A,T,x1lim,x2lim);
plot(xs(:,1),xs(:,2),'k--','LineWidth',5)

% Nullclines
x1node = nodeunif(100,x1lim(1),x1lim(2));
x2null = b0/b2 - (b1/b2)*x1node;
plot(x1node,x2null,'g','LineWidth',2)
plot([a1 a1],x2lim,'b','LineWidth',2)
plot(x1lim,[a2 a2],'b','LineWidth',2)

% Velocity Field
odefield(f,x1lim,x2lim)


%% SOLVE ODE USING ODERK4

% Set time discretization
N  = 100;           % number of time nodes

% Solve for different initial values and plot solution
for iv=1:5
  % Set initial values and horizon
  switch iv
    case 1          % initial value left of separatrix
      T  = 1;                         % time horizon
      x0 = [0; 2];                    % initial values
    case 2          % initial value left of separatrix
      T  = 3;                         % time horizon
      x0 = [0.9; 4.5];                  % initial values
    case 3          % initial value near unstable steady-state in direction of unstable eigenvector
      T  = 6;                         % time horizon
      J = [A(2)-a2 A(1)-a1;-b1 -b2];
      [V,D] = eig(J);
      x0 = A - 0.01*V(:,2);         % initial values
    case 4          % initial value on separatrix
      T  = 5;                         % time horizon
      x0 = [1; 1.5];                  % initial values
    case 5          % initial value right of separatrix
      T  = 5;                         % time horizon
      x0 = [2; 1];                    % initial values
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
text(A(1)+0.10,A(2)+0.3,'A','FontSize',14)
bullet(B(1),B(2),20)
text(B(1)+0.01,B(2)+0.3,'B','FontSize',14)


% %% Save Plots as EPS Files
% printfigures(mfilename,4)