% ICSOLVE Solves continuous time impulse control models
% USAGE
%   [cv,basis,x]=icsolve(model,x,n,type);
% INPUTS
%   model   : a structure variable (defined below)
%   x       : initial control values (2 by 2)
%   n       : degree of approximation
%   type    : the type of basis functions to use (default='cheb')
% OUTPUTS
%   cv      : solution coefficients (n by 1)
%   basis  : function definition structure
%   x       : optimal control values (2 by 2)
%
% The model structure should contain the following fields:
%   func    : model function file name (see below)  
%   params  : additional parameters to pass to model function
%   xindex  : defines the nature of the control values, see below (2 by 2)
%   R       : marginal rewards (2 by 1)
%   F       : fixed costs
%
%   The problem is defined on [a,b]. 
%   If a is a trigger, the process jumps to A when a is hit
%   If b is a trigger, the process jumps to B when b is hit
%   where a<=A<=B<=b
% x is a 2 by 2 matrix composed of [a A;b B]
%
% The 2 by 2 variable xindex defines the nature of the control x
%   i j xindex(i,j)  Condition                      
%   1 1    0         a is not a trigger
%   1 2    0
%                    a is a trigger and
%   1 1    1            a is a not a choice variable
%   1 1    2            a is a choice variable
%   1 2    0            A is a not a choice variable
%   1 2    1            A is a choice variable
%
%   2 1    0         b is not a trigger
%   2 2    0
%                    b is a trigger and
%   2 1    1            b is a not a choice variable
%   2 1    2            b is a choice variable
%   2 2    0            B is a not a choice variable
%   2 2    1            B is a choice variable
%
% The model function file should have the format
%    function out=func(flag,s,additional parameters)
%      switch flag
%        case 'f'
%          Return a matrix f representing the reward function
%        case 'mu'
%          Return a matrix mu representing the drift function
%              of the state transition equation
%        case 'sigma'
%          Return a matrix sigma representing the diffusion function
%              of the state transition equation
%        case 'rho'
%          Return a vector representing the (state contingent)
%              discount rates
%        case 'R+'
%          Return the reward associated with trigger s(1) and target s(2) (s(1)<s(2))
%            [R(s(1),s(2)) R_S0(s(1),s(2)) R_S1(s(1),s(2)) R_S0S0(s(1),s(2))+R_S0S1(s(1),s(2))]
%        case 'R-'
%          Return the reward associated with trigger s(1) and target s(2) (s(1)>s(2))
%            [R(s(1),s(2)) R_S0(s(1),s(2)) R_S1(s(1),s(2)) R_S0S0(s(1),s(2))+R_S0S1(s(1),s(2))]
%      end
% Each of these should be n by 1 where n=# of rows in s

% Copyright (c) 1997-2002, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function [cv,basis,x]=icsolve(model,x,n,type)

  if ~exist('type','var'), type='cheb'; end

% Unpack model variables
  func=model.func;
  params=model.params;
  xindex=model.xindex;
  F=model.F;

  % Get basis matrices and nodal points
  m=sum(xindex(:,1)>0);
  basis=fundefn(type,n-m,0,1);
  z=funnode(basis);
  basis=fundefn(type,n,0,1);
  Phi0=funbase(basis,z,0);
  Phi1=funbase(basis,z,1);
  Phi2=funbase(basis,z,2);

  phi0=funbase(basis,[0;1],0);
  phi1=funbase(basis,[0;1],1);
  phi2=funbase(basis,[0;1],2);
 
  % Define xindex variable
  if F(1)==0; xindex(1,2)=0; end
  if F(2)==0; xindex(2,2)=0; end
  xindex(:,1)=xindex(:,1)-1;

  % Call root finding algorithm (broyden) and get value function coefficients
  y=x(xindex>0);
  y=broyden('icres',y,x,func,params,xindex,basis,z,Phi0,Phi1,Phi2,phi0,phi1,phi2,F);
  x(xindex>0)=y;  x(F==0,2)=x(F==0,1);
  [e,cv]=icres(y,x,func,params,xindex,basis,z,Phi0,Phi1,Phi2,phi0,phi1,phi2,F);

  % Adjust basis 
  basis0=basis;
  basis=fundefn(type,n,x(1,1),x(2,1));