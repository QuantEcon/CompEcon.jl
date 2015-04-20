function demslv07

%% DEMSLV07 Linear complementarity problem methods

% Preliminary tasks
demosetup(mfilename)

% Generate test problem data
n  = 8;
z = randn(n,2)-1;
a = min(z,[],2);
b = max(z,[],2);
q  = randn(n,1);
M  = randn(n,n);
M  = -M'*M;
xinit = randn(n,1);

% Print table header
fprintf('Hundreds of seconds required to solve randomly generated linear\n')
fprintf('complementarity problem on R^8 using Newton and Lemke methods\n')
fprintf('Algorithm           Time      Norm\n');

% Solve by applying Newton method to minmax formulation
tic
optget('lcpsolve','type','minmax');
[x,z] = lcpsolve(M,q,a,b,xinit);
fprintf('Newton minmax     %6.2f  %8.0e\n',100*toc,norm(minmax(x,a,b,z)))

% Solve by applying Newton method to semismooth formulation
tic
optget('lcpsolve','type','smooth');
[x,z] = lcpsolve(M,q,a,b,xinit);
fprintf('Newton semismooth %6.2f  %8.0e\n',100*toc,norm(minmax(x,a,b,z)))

% Solve using Lemke's algorithm
tic
[x,z] = lcplemke(M,q,a,b,xinit);
fprintf('Lemke             %6.2f  %8.0e\n',100*toc,norm(minmax(x,a,b,z)))