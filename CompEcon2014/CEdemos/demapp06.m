function demapp06

%% DEMAPP06 Chebychev and cubic spline derivative approximation errors

% Preliminary tasks
demosetup(mfilename)

% Set degree of approximation and endpoints of approximation interval
a =  -1;                            % left endpoint
b =   1;                            % right endpoint
n =  10;                            % order of interpolatioin

% Construct refined uniform grid for error ploting
x = nodeunif(1001,a,b);

% Compute actual and fitted values on grid
[y,d,s] = f(x);                     % actual

% Construct and evaluate Chebychev interpolant
chebbasis = fundefn('cheb',n,a,b);  % chose basis ftions
c     = funfitf(chebbasis,@f);      % compute basis coefficients
ycheb = funeval(c,chebbasis,x);     % values
dcheb = funeval(c,chebbasis,x,1);   % first derivative
scheb = funeval(c,chebbasis,x,2);   % second derivative

% Construct and evaluate cubic spline interpolant
splibasis = fundefn('spli',n,a,b);  % chose basis ftions
c     = funfitf(splibasis,@f);      % compute basis coefficients
ycspl = funeval(c,splibasis,x);     % values
dcspl = funeval(c,splibasis,x,1);   % first derivative
scspl = funeval(c,splibasis,x,2);   % second derivative

% Plot function approximation error
figure
subplot(2,1,1); plot(x,y-ycheb,'r'); ylabel('Chebychev');
subplot(2,1,2); plot(x,y-ycspl,'m'); ylabel('Cubic Spline');
subplot(2,1,1); title('Function Approximation Error');
subplot(2,1,2); xlabel('x');

% Plot first derivative approximation error
figure
subplot(2,1,1); plot(x,d-dcheb,'r'); ylabel('Chebychev');
subplot(2,1,2); plot(x,d-dcspl,'m'); ylabel('Cubic Spline');
subplot(2,1,1); title('First Derivative Approximation Error');
subplot(2,1,2); xlabel('x');

% Plot second derivative approximation error
figure
subplot(2,1,1); plot(x,s-scheb,'r'); ylabel('Chebychev');
subplot(2,1,2); plot(x,s-scspl,'m'); ylabel('Cubic Spline');
subplot(2,1,1); title('Second Derivative Approximation Error');
subplot(2,1,2); xlabel('x');

% Save Plots as EPS Files
printfigures(mfilename,3)


%% Function to be approximated
function [y,d,s] = f(x)
y = exp(-x);
d = -exp(-x);
s = exp(-x);