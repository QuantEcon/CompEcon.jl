function demslv14

%% DEMSLV14 Spacial Equilibrum Model

% Preliminary tasks
demosetup(mfilename)

% Model parameters
as = [9; 3; 18];
bs = [1; 2; 1];
ad = [42; 54; 51];
bd = [3; 2; 1];
c = [0 3 9;3 0 3;6 3 0];
a = zeros(9,1);
b = inf*ones(9,1);

% Autarkic solution
% x = (ad-as)./(bs+bd);
% x = diag(x);
% x = reshape(x,9,1);

% Solve spatial equilibrium model
optset('ncpsolve','ShowIters',1)
x = zeros(9,1);
[x,fval] = ncpsolve(@f,a,b,x,as,bs,ad,bd,c);

% Output
norm(fval)
x = reshape(x,3,3)
p = as + bs.*sum(x,2)


function [fval,fjac] = f(x,as,bs,ad,bd,c)
x = reshape(x,3,3);
ps = as + bs.*(sum(x,2));
pd = ad - bd.*(sum(x,1))';
ps = ps(:,ones(1,3));
pd = pd(:,ones(1,3))';
fval = pd - ps - c;
fval = reshape(fval,9,1);
fjac = [];