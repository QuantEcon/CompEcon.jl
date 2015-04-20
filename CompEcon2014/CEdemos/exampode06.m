% exampode06

% Script file for predator-prey ODE example in lecture notes.

close all
clear all

% Model Parameters
a1 = 0.40;      % natural death rate of predators
a2 = 0.01;      % growth rate of predators per predation
b1 = 1.00;      % natural growth rate of prey
b2 = 0.02;      % death rate of prey per predation
h  = 0;         % human hunting of predators

% Velocity Function
f = @(x) [-a1*x(1,:)+a2*x(1,:).*x(2,:)-h;b1*x(2,:)-b2*x(1,:).*x(2,:)];