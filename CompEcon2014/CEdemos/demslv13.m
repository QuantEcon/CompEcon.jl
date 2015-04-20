function demslv13


%% DEMSLV13 Simple nonlinear complementarity problem

% Preliminary tasks
demosetup(mfilename)

a = [0;0]; b = [1;1];
x = [.5; .5];
optset('ncpsolve','ShowIters',1)
[x,f] = ncpsolve(@func,a,b,x)


function [fval,fjac] = func(x)
fval = [1+x(1)*(x(2)-2*x(1)^2-1); 2*x(1)^2-x(2)];
fjac = [];
% fjac = [-6*x(1)^2+x(2)-1 x(1); 4*x(1) -1];