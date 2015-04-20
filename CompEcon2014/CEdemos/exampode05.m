% exampode05

% Script file for competitive storage ODE example in lecture notes.

close all
clear all

% Model Parameters
r   = 0.1;                              % interest rate
k   = 0.5;                              % unit cost of storage
eta = 5;                                % demand elasticity

% Velocity Function
f = @(x) [-x(2,:).^(-eta);r*x(2,:)+k];