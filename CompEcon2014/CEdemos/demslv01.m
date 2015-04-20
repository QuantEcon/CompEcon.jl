function demslv01

%% DEMSLV01 Compute root of f(x)=exp(-x)-1
%
% Compute root of f(x)=exp(-x)-1 using Newton and secant methods. Initial 
% value generated randomly. True root is x=0.  

% Preliminary tasks
demosetup(mfilename)

% Randomly generate starting point
xinit = 10*randn;

% Print table header
fprintf('Hundreds of seconds required to compute root of exp(-x)-1,\n')
fprintf('via Newton and Broyden methods, starting at x=%4.2f.\n',xinit)
fprintf('Method      Time   Norm of f   Final x\n')

% Compute root using Newton method
tic;
x = newton(@f,xinit);
fprintf('Newton  %8.2f    %8.0e     %5.2f\n',100*toc,norm(f(x)),x)

% Compute root using Broyden method
tic;
x = broyden(@f,xinit);
fprintf('Broyden %8.2f    %8.0e     %5.2f\n',100*toc,norm(f(x)),x)


%% Function
function [fval,fjac] = f(x)
fval = exp(-x)-1;
fjac = -exp(-x);