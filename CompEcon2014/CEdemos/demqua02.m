function demqua02

%% DEMQUA02 Compute expectation of function of bivariate normal vector using Monte
% Carlo and Gaussian quadrature

% Preliminary tasks
demosetup(mfilename)

mu    = [0; 0];                  % enter mean vector
sigma = [1 0.5; 0.5 1];          % enter variance matrix

% Monte Carlo integration
n = 50000;                       % enter number of replicates
x = montnorm(n,mu,sigma);        % generate pseudo-random sequence
yexp = sum(f(x))/n;              % perform monte carlo integration of f
fprintf('Monte Carlo Integration:  %10.3f\n',yexp)

% Gaussian integration
n = [21 21];                     % enter order of approximation
[x,w] = qnwnorm(n,mu,sigma);     % compute Gaussian normal nodes and weights
yexp = w'*f(x);                  % perform Gaussian integration of f
fprintf('Guassian Quadrature:      %10.3f\n',yexp)

function y=f(x)
y = exp(-x(:,1)).*cos(x(:,2).^2);