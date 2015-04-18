% BAW Barone-Adesi/Whaley American option pricing model
% USAGE
%   [V,v,Sstar]=baw(sigma,S,K,r,delta,T,put);
% INPUTS
%   sigma : volatility
%   S     : price of underlying asset
%   K     : strike price
%   r     : interest rate on a discount bond
%   delta : dividend rate
%   T     : time to expiration
%   put   : 1 for puts, 0 for calls
% OUTPUTS
%   V     : American premium
%   v     : European premium (Black-Scholes)
%   Sstar : approximate early exercise boundary
%  
% The approximation can be written as
%   V(S,T)=v(S,T) + A*S^beta
% where beta solves
%   beta*(beta-1)*sigma^2/2 + (r-delta)*beta - r/(1-exp(-rT)) = 0
% and
% A and S* are determined by the boundary conditions
%   V(S*,T)=S-K, V_S(S*,T)=1 for calls
%   V(S*,T)=K-S, V_S(S*,T)=-1 for puts

% Copyright (c) 1997-2002, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function [V,v,Sstar]=baw(sigma,S,K,r,delta,T,put);

  sigma(sigma<=0|sigma>1e8)=NaN;

  tol=K*5e-5;
  call=(put<=0)-(put>0);          % call=1, put=-1
  B=exp(-r.*T);
  D=exp(-delta.*T);
  sigmat=call.*sigma.*sqrt(T);

  sigma2=sigma.*sigma;
  temp1=2*(r-delta)./sigma2;
  temp2=2*r./sigma2./(1-B);
  beta=(1-temp1+call.*sqrt((1-temp1).^2+4*temp2))/2;

  % Determine the early exercise price using Newton's method
  kk=K./(1-1./beta);
  factor=(r-delta).*T./sigmat+0.5*sigmat;
  factor1=B./2.5066282746310001./sigmat.*kk;  % sqrt(2*pi)=2.5066282746310001
  factor2=D./2.5066282746310001./sigmat;
  S0=0;
  Sstar=kk;
  while any(abs(Sstar-S0)>tol)
    S0=Sstar;
    d1=log(Sstar./K)./sigmat+factor;
    d2=d1-sigmat;
    temp1=1-D.*cdfn(d1);
    temp2=1-B.*cdfn(d2);
    res=Sstar.*temp1-kk.*temp2;
    dres=temp1+factor1./Sstar.*exp(-0.5*d2.*d2)-factor2.*exp(-0.5*d1.*d1);
    Sstar=Sstar-res./dres;
  end
  
  % Compute the option value
  A=call.*temp1.*(Sstar.^(1-beta))./beta;
  v=bs(sigma,S,K,r,delta,T,put);
  V=v+A.*(S.^beta);
  
  % Adjust for early exercise
  indc=S>=Sstar & call==1;
  indp=S<=Sstar & call==-1;
  V=indc.*(S-K) + indp.*(K-S) + ~(indp+indc).*V;
