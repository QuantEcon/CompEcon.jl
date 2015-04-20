function demapp02

%% DEMAPP02 Approximating functions on R^2

% This m-file illustrates how to use CompEcon Toolbox routines to construct
% and operate with an approximant for a function defined on a rectangle in 
% R^2.

% In particular, we construct an approximant for f(x1,x2)=cos(x1)/exp(x2)
% on [-1,1]X[-1,1].  The function used in this illustration posseses a 
% closed-form, which will allow us to measure approximation error precisely.
% Of course, in practical applications, the function to be approximated 
% will not possess a known closed-form.

% In order to carry out the exercise, one must first code the function to 
% be approximated at arbitrary points.  The required code is presented at
% the end of this m-file (see below). Let's begin:

% Preliminary tasks
demosetup(mfilename)

% Set the endpoints of approximation interval:
a = [-1 -1];                          % left endpoints
b = [ 1  1];                          % right endpoints

% Choose an approximation scheme. In this case, let us use an 11 by 11 
% Chebychev approximation scheme:
n = [11 11];                          % order of approximation
basis = fundefn('cheb',n,a,b);        % define basis

% Compute the basis coefficients c.  There are various way to do this:
% One may use funfitf:
c = funfitf(basis,@f);  

% ... or one may compute the standard approximation nodes x and corresponding
% function values y and use funfitxy:
x = funnode(basis);
y = f(x);
c = funfitxy(basis,x,y);

% ... or one compute the standard approximation nodes x, corresponding
% function values y, and the interpolation matrix phi, and solve the 
% interpolation equation directly using the backslash operator:
x = funnode(basis);
y = f(x);
phi = funbase(basis);
c = phi\y;

% Having computed the basis coefficients, one may now evaluate the 
% approximant at any point x using funeval:
x = [0 0];
y = funeval(c,basis,x);
fprintf('The exact and approximate value of f at x=[0 0] are\n')
fprintf('%4.0f  %20.15f\n\n',0,y)

% ... one may also evaluate the approximant's first partial derivatives 
% at x:
d1 = funeval(c,basis,x,[1 0]);
d2 = funeval(c,basis,x,[0 1]);
fprintf('The exact and approximate partial derivatives of f w.r.t. x1 at x=[0 0] are\n');
fprintf('%4.0f  %20.15f\n\n',1,d1)
fprintf('The exact and approximate partial derivatives of f w.r.t. x2 at x=[0 0] are\n');
fprintf('%4.0f  %20.15f\n\n',0,d2)

% ... one may also evaluate the approximant's second own partial and cross 
% partial derivatives at x:
d11 = funeval(c,basis,x,[2 0]);
d22 = funeval(c,basis,x,[0 2]);
d12 = funeval(c,basis,x,[1 1]);
fprintf('The exact and approximate second partial derivatives of f w.r.t. x1 at x=[0 0] is\n');
fprintf('%4.0f  %20.15f\n\n',0,d11)
fprintf('The exact and approximate second partial derivatives of f w.r.t. x2 at x=[0 0] is\n');
fprintf('%4.0f  %20.15f\n\n',0,d22)
fprintf('The exact and approximate second cross partial derivatives of f at x=[0 0] is\n');
fprintf('%4.0f  %20.15f\n\n',-1,d12)

% One may evaluate the accuracy of the Chebychev polynomial approximant by 
% computing the approximation error on a highly refined grid of points:
nplot = [101 101];                    % chose grid discretization
[xgrid,xcoord] = nodeunif(nplot,a,b); % generate refined grid for plotting
yapp = funeval(c,basis,xgrid);        % approximant values at grid nodes
yact = f(xgrid);                      % actual function values at grid points
error = reshape(yapp-yact,nplot(1),nplot(2));
figure; surf(xcoord{1},xcoord{2},error,'FaceColor','interp','EdgeColor','interp');
xlabel('x1'); ylabel('x2'); zlabel('error');
title('Chebychev Approximation Error')
% The plot indicates that an order 11 by 11 Chebychev approximation scheme 
% produces approximation errors no bigger in magnitude than 1x10^-10. 

% Let us repeat the approximation exercise, this time constructing an order
% 21 by 21 cubic spline approximation scheme:
n = [21 21];                          % order of approximation
basis = fundefn('spli',n,a,b);        % define basis
c = funfitf(basis,@f);                % compute basis coefficients
yapp = funeval(c,basis,xgrid);        % approximant values at grid nodes
error = reshape(yapp-yact,nplot(1),nplot(2));
figure; surf(xcoord{1},xcoord{2},error,'FaceColor','interp','EdgeColor','interp');
xlabel('x1'); ylabel('x2'); zlabel('error');
title('Cubic Spline Approximation Error')
% The plot indicates that an order 21 by 21 cubic spline approximation
% scheme produces approximation errors no bigger in magnitude than 1x10^-6.

% Save Plots as EPS Files
printfigures(mfilename,2)


%% Function to be approximated
function y = f(x)
y = sin(x(:,1))./exp(x(:,2));