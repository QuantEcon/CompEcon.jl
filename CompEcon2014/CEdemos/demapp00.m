function demapp00

%% DEMAPP00 Approximating Using the CompEcon Toolbox

% Preliminary tasks
demosetup(mfilename)


%% Univariate Interpolation

% Define Function and Derivative
f1 = @(x) exp(-2*x);
d1 = @(x) -2*exp(-2*x);

% Fit approximant
n = 10; a = -1; b = 1;
basis = fundefn('cheb',n,a,b);
c = funfitf(basis,f1);

% Graph approximation error for function and derivative
x = nodeunif(1001,a,b);
yact = f1(x);
dact = d1(x);
yfit = funeval(c,basis,x);
dfit = funeval(c,basis,x,1);

% Nice plot of function approximation error
figure
hold on
plot(x,yfit-yact);
plot(x,0*x,'k--','linewidth',2);
set(gca,'xtick',[-1 0 1])
xlabel('x')
ylabel('Error')

% Nice plot of derivative approximation error
figure
hold on
plot(x,dfit-dact);
plot(x,0*x,'k--','linewidth',2);
set(gca,'xtick',[-1 0 1])
xlabel('x')
ylabel('Error')


%% Bivariate Interpolation

% Define Function
f2 = @(x) cos(x(:,1))./exp(x(:,2));

% Set degree and domain of interpolation
n = [7 7];
a = [ 0 0];
b = [ 1 1];
basis = fundefn('cheb',n,a,b);
c = funfitf(basis,f2);

% Nice plot of function approximation error
figure
nplot = [101 101];
[x,xcoord] = nodeunif(nplot,a,b);
yfit = funeval(c,basis,x);
error = reshape(yfit-f2(x),nplot);
surf(xcoord{1},xcoord{2},error,'FaceColor','interp','EdgeColor','interp');
xlabel('x_1'); ylabel('x_2'); zlabel('Error');
title('Chebychev Approximation Error')

% Compute partial Derivatives
x = [0.5 0.5]
f1  = funeval(c,basis,x,[1 0])
f2  = funeval(c,basis,x,[0 1])
f11 = funeval(c,basis,x,[2 0])
f12 = funeval(c,basis,x,[1 1])
f22 = funeval(c,basis,x,[0 2])


%% Save Plots as EPS Files
printfigures(mfilename,3)