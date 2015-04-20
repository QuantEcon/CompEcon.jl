function demslv04

%% DEMSLV04 Compute fixedpoint of f(x1,x2)= [x1^2+x2^3;x1*x2-0.5]
%
% Compute fixedpoint of f(x1,x2)= [x1^2+x2^3;x1*x2-0.5] using Newton, 
% Broyden, and function iteration methods. Initial values generated 
% randomly.  Some algorithms may fail to converge, depending on the
% initial value.  True fixedpoint is x1=-0.09 x2=-0.46.  

% Preliminary tasks
demosetup(mfilename)

% Randomly generate starting point
xinit = randn(2,1);

% Print table header
fprintf('Hundreds of seconds required to compute fixed-point of g(x1,x2)=[x1^2+x2^3;x1*x2-0.5]\n')
fprintf('using Newton, Broyden, and function iteration methods, starting at x1=%4.2f x2=%4.2f\n',xinit)
fprintf('Method      Time   Norm of f        x1     x2\n')

% Compute fixed-point using Newton method
optset('newton','maxit',1500);
tic;
x = newton(@f,xinit);
fprintf('Newton  %8.2f    %8.0e     %5.2f  %5.2f\n',100*toc,norm(f(x)),x)

% Compute fixed-point using Broyden method
optset('broyden','maxit',1500);
tic;
x = broyden(@f,xinit);
fprintf('Broyden %8.2f    %8.0e     %5.2f  %5.2f\n',100*toc,norm(f(x)),x)

% Compute fixed-point using function iteration
optset('fixpoint','maxit',1500);
tic;
x = fixpoint(@g,xinit);
fprintf('Function%8.2f    %8.0e     %5.2f  %5.2f\n',100*toc,norm(f(x)),x)


%% Function
function gval = g(x)
gval = [x(1)^2 + x(2)^3; x(1)*x(2) - 0.5];


%% Equivalent Rootfinding Formulation
function [fval,fjac] = f(x)
fval = [x(1) - x(1)^2 - x(2)^3; x(2) - x(1)*x(2) + 0.5];
fjac = [1-2*x(1) -3*x(2)^2; -x(2) 1-x(1)];