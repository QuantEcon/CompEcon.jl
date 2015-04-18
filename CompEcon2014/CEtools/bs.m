%% BS 
%
%  Computes the value of a European put or call option using the 
%  Black-Scholes option pricing formula
%
%  Usage
%    V = bs(sigma,S,K,r,delta,T,put)
%  Input
%    sigma : annaulized volatility
%    S     : price of underlying asset
%    K     : strike price
%    r     : annualized risk-free rate
%    delta : annualized dividend rate
%    T     : time to expiration in years
%    put   : 1 if put option, 0 if call option
%  Output
%    V     : option premium

%  Copyright(c) 1997-2010
%   Mario J. Miranda - miranda.4@osu.edu
%   Paul L. Fackler  - paul_fackler@ncsu.edu

function [V,d1] = bs(sigma,S,K,r,delta,T,put)

sigma = sigma.*sqrt(T);
S  = exp(-delta.*T).*S;
K  = exp(-r.*T).*K;
d1 = log(S./K)./sigma + sigma/2;
V  = S.*cdfn(d1)-K.*cdfn(d1-sigma);
V  = (~put).*V + put.*(V-S+K);