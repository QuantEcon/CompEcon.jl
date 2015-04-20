% exampode07

% Script file for fisheries model ODE example in lecture notes.

close all
clear all

% Model Parameters
alpha = 0.5;
beta  = 2.75;
phi   = 0.05;
delta = 10;

% Velocity Function
s = @(x) x(1,:);
k = @(x) x(2,:);
y = @(x) s(x)./(alpha+beta*s(x).*k(x));
p = @(x) alpha./(alpha+beta*s(x).*k(x));
f = @(x) [(1-s(x)).*s(x)-k(x).*y(x); delta*(0.5*p(x).*y(x)-phi)];