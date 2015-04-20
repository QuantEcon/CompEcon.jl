function demapp01

%% DEMAPP01 Approximating functions on R

% This m-file illustrates how to use CompEcon Toolbox routines to construct
% and operate with an approximant for a function defined on an interval of 
% the real line.

% In particular, we construct an approximant for f(x)=exp(-x) on the 
% interval [-1,1].  The function used in this illustration posseses a 
% closed-form, which will allow us to measure approximation error precisely.
% Of course, in practical applications, the function to be approximated 
% will not possess a known closed-form.

% In order to carry out the exercise, one must first code the function to 
% be approximated at arbitrary points.  The required code is presented at
% the end of this m-file (see below). Let's begin:

% Preliminary tasks
demosetup(mfilename)

% Set the endpoints of approximation interval:
a =  -1;                            % left endpoint
b =   1;                            % right endpoint

% Choose an approximation scheme. In this case, let us use an order 10 
% Chebychev approximation scheme:
n = 10;                             % order of approximation
basis = fundefn('cheb',n,a,b);      % define basis

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
x = 0;
y = funeval(c,basis,x);
disp('The approximate value of exp(-x) at x=0 is'); disp(y);
disp('The ''exact'' value of exp(-x) at x=0 is'); disp(1);

% ... one may also evaluate the approximant's first and second derivatives 
% at x:
d1 = funeval(c,basis,x,1);
d2 = funeval(c,basis,x,2);
disp('The approximate first derivative of exp(-x) at x=0 is'); disp(d1);
disp('The ''exact'' first derivative of exp(-x) at x=0 is'); disp(-1);
disp('The approximate second derivative of exp(-x) at x=0 is'); disp(d2);
disp('The ''exact'' second derivative of exp(-x) at x=0 is'); disp(1);

% ... and one may even evaluate the approximant's definite integral between 
% the left endpoint a and x:
int = funeval(c,basis,x,-1);
disp('The approximate integral of exp(-x) between x=-1 and x=0 is'); disp(int);
disp('The ''exact'' integral of exp(-x) between x=-1 and x=0 is'); disp(exp(1)-1);

% One may evaluate the accuracy of the Chebychev polynomial approximant by 
% computing the approximation error on a highly refined grid of points:
ngrid = 5001;                       % number of grid nodes
xgrid = nodeunif(ngrid,a,b);        % generate refined grid for plotting
yapp = funeval(c,basis,xgrid);      % approximant values at grid nodes
yact = f(xgrid);                    % actual function values at grid points
figure                              % plot approximation error
hold on
plot(xgrid,yapp-yact)               
plot(xgrid,zeros(ngrid,1),'k--','LineWidth',2) 
xlabel('x'); ylabel('Error'); 
title('Chebychev Approximation Error for exp(-x)\newline')
% The plot indicates that an order 10 Chebychev approximation scheme, 
% produces approximation errors no bigger in magnitude than 6x10^-10. The 
% approximation error exhibits the "Chebychev equioscillation property", 
% oscilating relatively uniformly throughout the approximation domain.  
% This commonly occurs when function being approximated is very smooth, as
% is the case here; but should not be expected when the function is not 
% smooth.  Further notice how the approximation error is exactly 0 at the 
% approximation nodes --- which is true by contruction.

% Let us repeat the approximation exercise, this time constructing a 
% 21-function cubic spline approximant:
n = 21;                             % order of approximation
basis = fundefn('spli',n,a,b);      % define basis
c = funfitf(basis,@f);              % compute basis coefficients
yapp = funeval(c,basis,xgrid);      % approximant values at grid nodes
figure                              % plot approximation error
hold on
plot(xgrid,yapp-yact)              
plot(xgrid,zeros(ngrid,1),'k--','LineWidth',2) 
xlabel('x'); ylabel('Error');
title('Cubic Spline Approximation Error for exp(-x)')
% The plot indicates that an order 21 cubic spline approximation scheme
% produces approximation errors no bigger in magnitude than 1.2x10^-6, about 
% four orders of magnitude worse than with Chebychev polynomials.

% Let us repeat the approximation exercise, this time constructing a 
% 31-function linear spline approximant:
n = 31;                             % order of approximation
basis = fundefn('spli',n,a,b,1);    % define basis functions
c = funfitf(basis,@f);              % compute basis coefficients
yapp = funeval(c,basis,xgrid);      % fitted values at grid nodes
figure                              % plot approximation error
hold on
plot(xgrid,yapp-yact)              
plot(xgrid,zeros(ngrid,1),'k--','LineWidth',2) 
xlabel('x'); ylabel('Error');
title('Linear Spline Approximation Error for exp(-x)')
% The plot indicates that an order 41 linear spline approximation scheme
% produces approximation errors no bigger in magnitude than 1.5x10^-2, about 
% five orders of magnitude worse than with Chebychev polynomials.

% Save Plots as EPS Files
% printfigures(mfilename,3)


%% Function to be approximated
function y = f(x)
y = exp(-x);