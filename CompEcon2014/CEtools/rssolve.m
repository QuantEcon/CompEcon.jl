% RSSOLVE Solves continuous time regime switching models
% Solves m-regime switching model with one continuous state variable
% USAGE
%   [cv,basis,x]=rssolve(model,x,n,type);
% INPUTS
%   model          : a structure variable (defined below)
%   x              : initial control values (px1)
%   n              : degree of approximation (m-vector)
%   type           : the type of basis functions to use (default='cheb')
% OUTPUTS
%   cv             : 1 by m cell array of solution coefficients 
%                      with cv{i} n(i) by 1
%   basis         : 1 by m cell array of function definition structures
%   x              : optimal switch points (p by 1)
%
% The model structure should contain the following fields:
%   func     : model function file (see below)  
%   params   : additional parameters to pass to model function file
%   xindex   : defines the nature of the control values
%                (p by 5) see below 
%
%   Regime i has a no switch region [x(1,i),x(2,i)]. 
%   When S hits x(j,i) the discrete regime switches to xindex(j,i)
%   If x(j,i) is a boundary for S x(j,i)=i and no switch occurs
%   
%   xindex defines the topology of switches.
%   If there is a point at which i switches to j, set
%     xindex(k,1)=i and xindex(k,2)=j 
%   and define the type of side constraints to hold at that point:
%       xindex(k,3)=>0 implies   V(S,i) =   V(S,j)+reward(k,1)
%       xindex(k,4)=>0 implies  V'(S,i) =  V'(S,j)+reward(k,2)
%       xindex(k,5)=>0 implies V''(S,i) = V''(S,j)+reward(k,3)
%   where reward is returned by the model function file.
%   To impose a side condition on V(S,i) alone, set
%     xindex(k,1)=i and xindex(k,2)=0 
%   and define the type of side constraints to hold at that point:
%       xindex(k,3)=>0 implies   V(S,i) = reward(k,1)
%       xindex(k,4)=>0 implies  V'(S,i) = reward(k,2)
%       xindex(k,5)=>0 implies V''(S,i) = reward(k,3)
%
% The model function file should have the format
%    out1=func(flag,s,additional parameters)
%      switch flag
%        case 'f'
%          Return a matrix f representing the reward function
%        case 'g'
%          Return a matrix g representing the drift function
%              of the state transition equation
%        case 'sigma'
%          Return a matrix sigma representing the diffusion function
%              of the state transition equation
%        case 'rho'
%          Return a vector representing the (state contingent)
%              discount rates
%        case 'reward'
%          Return a k x 3 matrix of jump rewards and derivatives
%              when passed a k x 1 matrix of state values
%      end
%
% USER OPTIONS (SET WITH OPTSET)
%   maxiters      : maximum number of iterations
%   tol           : convergence tolerance
%   showiters     : 0/1, 1 to display iteration results

% Copyright (c) 1997-2002,  Paul L. Fackler & Mario J. Miranda 
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function [cv,basis,x]=rssolve(model,x,n,type)

% Set default options
  maxiters      = optget('rssolve','maxiters',100);
  tol           = optget('rssolve','tol',sqrt(eps));
  showiters     = optget('rssolve','showiters',1);
  
  optset('broyden','maxiters',maxiters)
  optset('broyden','tol',tol)
  optset('broyden','showiters',showiters)

  if ~exist('type','var'), type='cheb'; end

% Unpack model variables
  func=model.func;
  params=model.params;
  xindex=model.xindex;

  if any(sum(xindex(:,3:5)==1,2)>2 |...
         sum(xindex(:,3:5)==2,2)>1)
    error('Incorrect specification of xindex');
  end

  m=max(max(xindex(:,[1 2])));
  if length(n)==1, n=n(ones(m,1)); end

  % Determine number of nodal points
  ncond=zeros(1,m);
  for k=1:size(xindex,1)
    i=xindex(k,1);
    j=xindex(k,2);
    s=sum(xindex(k,3:5)==1);
    if j==0
      ncond(i)=ncond(i)+s;
    else
      ncond(i)=ncond(i)+1;
      if s==2, ncond(j)=ncond(j)+1; end
    end
  end

  % Get basis matrices and nodal points
  for i=1:m
    basis{i}=fundefn(type,n(i)-ncond(i),0,1);
    z{i}=funnode(basis{i});
    basis{i}=fundefn(type,n(i),0,1);
    Phi0{i}=funbase(basis{i},z{i},0);
    Phi1{i}=funbase(basis{i},z{i},1);
    Phi2{i}=funbase(basis{i},z{i},2); 
    phi0{i}=funbase(basis{i},[0;1],0);
    phi1{i}=funbase(basis{i},[0;1],1);
    phi2{i}=funbase(basis{i},[0;1],2);
  end

  xchoice=logical(sum(xindex(:,3:5)==2,2));
  y=x(xchoice);
  y=broyden('rsres',y,x,func,params,xindex,xchoice,m,n,...
             basis,z,ncond,Phi0,Phi1,Phi2,phi0,phi1,phi2);
  x(xchoice)=y; 
  [e,c]=rsres(y,x,func,params,xindex,xchoice,m,n,...
             basis,z,ncond,Phi0,Phi1,Phi2,phi0,phi1,phi2);
 
  for i=1:m
    cv{i}=c(1:n(i));
    c(1:n(i))=[];
    j=xindex(:,1)==i |xindex(:,2)==i;
    a=min(x(j));
    b=max(x(j));
    basis{i}=fundefn(type,n(i),a,b);
  end
  
  optset('broyden','defaults')
 