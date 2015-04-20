% exampode03

% Script file for generic ODE example in lecture notes.

close all
clear all


% Model Parameters
A  = [1 12; -1 -6];     % velocity function parameters
b  = [-60; 36];         % velocity function parameters

% Velocity Function
f = @(x) [A(1,1)*x(1,:)+A(1,2)*x(2,:)+b(1); ...
          A(2,1)*x(1,:)+A(2,2)*x(2,:)+b(2)];
        
% Closed-Form Solution 
X = @(t) [12 - 48*exp(-2*t) + 42*exp(-3*t) ...
           4 + 12*exp(-2*t) - 14*exp(-3*t)];