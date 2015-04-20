% exampode12

% Script file for regioinal migration model ODE example in lecture notes.

close all
clear all


%% FORMULATION

% Model Parameters
theta = 0.5;        % inverse labor demand elasticity
k0    = 0.1;        % minimum migration cost
k1    = 0.1;        % migration cost parameter
k2    = 0.1;        % migration cost parameter
rho   = 0.05;       % discount rate
wbar  = 1;          % world wage rate

% Velocity Function
L = @(x) x(1,:);
B = @(x) x(2,:);
f = @(x) [L(x).*(sqrt(k1^2+4*k2*(B(x)-k0))-k1)/(2*k2); rho*B(x)-L(x).^(-theta)+wbar];