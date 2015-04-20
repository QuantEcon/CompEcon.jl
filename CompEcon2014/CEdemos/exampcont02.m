function exampcont02

% Script file for renewable resouce example in lecture notes.

close all
clear all


%% FORMULATION

% Model Parameters
alpha = 0.5;        % biological growth function scale factor
beta  = 0.5;        % biological growth function elasticity factor
kappa = 5;          % unit harvest cost scale factor
gamma = 1.5;        % unit harvest cost elasticity
eta   = 1.5;        % inverse elasticity of demand
rho   = 0.05;       % discount rate
sigma = 0.1;        % diffusion volatiity


%% SOLVE BELLMAN EQUATION USING SCSOLVE

% Begin typing here



%% Function File for scsolve
function out = func(flag,s,q,Vs,alpha,beta,kappa,gamma,eta,sigma)
k = kappa*s.^(-gamma);
switch flag
  case 'x'
    out = (Vs+k).^(-1/eta);
  case 'f'
    out = (1/(1-eta))*q.^(1-eta) - k.*q;
  case 'g'
    g = (alpha/beta)*s.*(1-s.^beta);
    out = g - q;
  case 'sigma'
    out = sigma*s;
end