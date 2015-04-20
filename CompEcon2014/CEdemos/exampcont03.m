function exampcont03

% Script file for economic growth example in lecture notes.

close all
clear all


%% FORMULATION

% Model Parameters
alpha =  0.4;       % capital share
delta =  0.05;      % capital depreciation rate
gamma =  1;         % technology mean reversion coefficient
theta =  2.0;       % relative risk aversion
rho   =  0.05;      % discount rate
sigma =  0.1 ;      % volatility


%% SOLVE BELLMAN EQUATION USING SCSOLVE

% Begin typing here




%% Function File for scsolve
function out = func(flag,s,q,Vs,alpha,delta,gamma,theta,sigma)
k = s(:,1);
y = s(:,2);
n = length(k);
switch flag
  case 'x'
    Vk = Vs(:,1);
    out = Vk.^(-1./theta);
  case 'f'
    out = (q.^(1-theta))./(1-theta);
  case 'g'
    out = [(y.*k.^alpha-delta*k-q)  gamma*(y-1)];
  case 'sigma'
    out = zeros(n,2,2);
    out(:,2,2) = sigma*sqrt(y);
end