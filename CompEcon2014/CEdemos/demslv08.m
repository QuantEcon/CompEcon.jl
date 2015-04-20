function demslv08

%% DEMSLV08 Nolinear complimentarty problem methods
%
% Solve nonlinear complementarity problem on R^2 using semismooth and
% minmax methods

% Preliminary tasks
demosetup(mfilename)
warning off

% Generate problem test data
z = randn(2,2)*2;
a = 1+min(z,[],2);
b = 1+max(z,[],2);
xinit = randn(2,1);

% Check Jacobian
[error,i,j] = checkjac(@f,xinit);

% Print table header
fprintf('Hundreds of seconds required to solve nonlinear complementarity\n')
fprintf('problem on R^2 using minmax and semismooth formulations, with\n')
fprintf('randomly generated bounds a = %4.2f %4.2f and b = %4.2f %4.2f\n',a,b)
fprintf('Algorithm           Time      Norm        x1     x2\n');

% Solve by applying Newton method to minmax formulation
tic
optset('ncpsolve','type','minmax');
optset('ncpsolve','maxit',1500);
[x,z] = ncpsolve(@f,a,b,xinit);
fprintf('Newton minmax     %6.2f  %8.0e     %5.2f  %5.2f\n',100*toc,norm(minmax(x,a,b,z)),x)

% Solve by applying Newton method to semismooth formulation
tic
optset('ncpsolve','type','ssmooth');
optset('ncpsolve','maxit',1500);
[x,z] = ncpsolve(@f,a,b,xinit);
fprintf('Newton semismooth %6.2f  %8.0e     %5.2f  %5.2f\n',100*toc,norm(minmax(x,a,b,z)),x)


%% Function file
function [fval,fjac] = f(x)
fval = [200*x(1)*(x(2)-x(1)^2)+1-x(1); 100*(x(1)^2-x(2))];
fjac = [200*(x(2)-x(1)^2)-400*x(1)^2-1 200*x(1); 200*x(1) -100];

% f = [1+x(1)*(x(2)-2*x(1)^2-1); 2*x(1)^2-x(2)];
% J = [-1+(x(2)-6*x(1)^2) x(1); 4*x(1) -1];