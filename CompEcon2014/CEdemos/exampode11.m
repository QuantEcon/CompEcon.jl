% exampode011

% Script file for Tobin's Q ODE example in lecture notes.

close all
clear all



%% FORMULATION

% Model Parameters
delta = 0.05;              % depreciation rate
b     = 5;                 % marginal adjustment cost
rho   = 0.05;              % interest rate

% Velocity Function
k = @(x) x(1,:);
q = @(x) x(2,:);
f = @(x) [k(x).*(max(q(x)-1,0)/b-delta); ...
  (delta+rho)*q(x)-1./k(x)-(q(x)-1).^2/(2*b)]% printfigures(mfilename,3)