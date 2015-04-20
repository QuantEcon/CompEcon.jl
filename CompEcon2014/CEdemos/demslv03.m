function demslv03

%% DEMSLV03 Compute fixedpoint of f(x) = x^0.5
%
% Compute fixedpoint of f(x) = x^0.5 using Newton, Broyden, and function 
% iteration methods. Initial values generated randomly. Some alrorithms 
% may fail to converge, depending on the initial value. True fixedpoint is 
% x=1.  

% Preliminary tasks
demosetup(mfilename)

% Randomly generate starting point
xinit = rand+0.5;

% Print table header
fprintf('Hundreds of seconds required to compute fixed-point of g(x)=sqrt(x)\n') 
fprintf('using Newton, Broyden, and function iteration methods, starting at\n')
fprintf('x=%4.2f\n',xinit)
fprintf('Method      Time   Norm of f         x\n')

% Compute fixed-point using Newton method
tic;
x = newton(@f,xinit);
fprintf('Newton  %8.2f    %8.0e     %5.2f\n',100*toc,norm(f(x)),x)

% Compute fixed-point using Broyden method
tic;
x = broyden(@f,xinit);
fprintf('Broyden %8.2f    %8.0e     %5.2f\n',100*toc,norm(f(x)),x)

% Compute fixed-point using function iteration
tic;
x = fixpoint(@g,xinit);
fprintf('Function%8.2f    %8.0e     %5.2f\n',100*toc,norm(f(x)),x)


%% Function
function gval = g(x)
gval = sqrt(x);


%% Equivalent Rootfinding Formulation
function [fval,fjac] = f(x)
fval = x-sqrt(x);
fjac = 1-0.5./sqrt(x);