% DEMAPP11
%
% Example of usage of fundefn and funeval with multiple dimension domain
% and range, with evaluation of multiple derivatives; also demonstrate use
% of funeval to evaluate integral


%% Differntiation Example
%
% f(x) = [x1.^2-x2 x1.*x2 x1-x2]
%
% dx = 2        dimension of function domain
% dy = 3        dimension of function range
% nx = 1000;    number of function/derivative evaluations
% n = [5 3];    number of basis functions per domain dimension
% p  = 4        number of function values/derivatives evaluated

% Approximation Structure
n = [5 3];
a = [0 0];
b = [1 1];
basis = fundefn('cheb',n,a,b);

% Define function
f = @(x) [x(:,1).^2-x(:,2) x(:,1).*x(:,2) x(:,1)-x(:,2)];

% Generate x and f(x) values and fit approximant
[x,xcoord] = nodeunif([200 80],a,b);
y = f(x);
c = funfitxy(basis,x,y);

% Compute approximate derivative values inputing x
y1 = funeval(c,basis,x,[1 1;0 1;1 0;2 0]);
size(y1)

% Compute approximate derivative values inputing xcoord
y2 = funeval(c,basis,xcoord,[1 1;0 1;1 0;2 0]);
norm(y2(:)-y1(:))

% Demonstrate use of funbasex
B = funbasex(basis,x,[0;1;2]);
y3 = funeval(c,basis,B,[1 1;0 1;1 0;2 0]);
norm(y3(:)-y1(:))


%% Integration Example
%
% Compute approximate integral of
%   f(x) = exp(x)
% on [0,1[

% Approximation Structure
n = 5;
a = 0;
b = 1;
basis = fundefn('cheb',n,a,b);

% Define function
f = @(x) exp(-x);

% Generate x and f(x) values and fit approximant
x = nodeunif(100,a,b);
y = f(x);
c = funfitxy(basis,x,y);

% Compute approximate integral
Y = funeval(c,basis,1,-1);
Ytrue = 1-1/exp(1);
norm(Y-Ytrue)