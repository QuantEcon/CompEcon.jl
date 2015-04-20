function demopt05

%% DEMOPT05 Optimization with qnewton

% Preliminary tasks
demosetup(mfilename)        

f = @(x) x^3-12*x^2+36*x+8;
x0 = 4;
x = qnewton(f,x0)
J = fdjac(f,x)
E = eig(fdhess(f,x))

f = @(x) 5-4*x(1)^2-2*x(2)^2-4*x(1)*x(2)-2*x(2);
x0 = [-1;1];
x  = qnewton(f,x0)
J = fdjac(f,x)
E = eig(fdhess(f,x))