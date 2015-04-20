% exampode09

% Script file for economic growth ODE example in lecture notes.

close all
clear all


%% FORMULATION

% Model Parameters
alpha = 0.5;               % capital share
delta = 0.05;              % depreciation rate
theta = 2.0;               % relative risk aversion
rho   = 0.05;              % discount rate

% Velocity Function
f = @(x) [x(1,:).^alpha-delta*x(1,:)-x(2,:); ...
         (alpha*x(1,:).^(alpha-1)-delta-rho).*x(2,:)/theta];