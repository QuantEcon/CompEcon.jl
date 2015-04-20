% exampode08 - nnonrenewable resource example in lecture notes.

close all
clear all

% Model Parameters
alpha = 2;          % biological growth function scale factor
kappa = 1;          % unit harvest cost scale factor
eta   = 2;          % inverse elasticity of demand
rho   = 0.05;       % discount rate

% Ancillary Functions
s    = @(x) x(1,:);                 % map x to s
q    = @(x) x(2,:);                 % map x to q
p    = @(x) q(x).^(-eta);           % inverse demand func
pder = @(x) -eta*q(x).^(-eta-1);    % inverse demand deriv
g    = @(s) alpha*s.*(1-s);         % biological growth func
gder = @(s) alpha*(1-2*s);          % biological growth deriv

% Velocity Function
f    = @(x) [g(s(x))-q(x); ...
   ((rho-gder(s(x))).*(p(x)-kappa))./pder(x)];