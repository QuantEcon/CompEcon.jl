function demopt08

%% DEMOPT08 Constrained optimization using fmincon

% Preliminary tasks
demosetup(mfilename)

x0 = [1;2];
[x,f,flag] = fmincon(@func,x0,[],[],[],[],[0;0],[],@cons)

function f = func(x)      
% f = -x(1)^2-x(2)^2; 
% f = -f;
f = -x(1)^2 - 2*x(2)^2 - 2*x(1)*x(2) + 1; 
f = -f; 
                                       
function [c,ceq] = cons(x)
% c = [];
% ceq = x(1)+x(2)-2;
ceq = [];
c = 1-x(1)-x(2);