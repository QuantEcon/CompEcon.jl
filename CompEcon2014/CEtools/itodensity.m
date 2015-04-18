% ITODENSITY Long-run densities for 1-D Ito processes
% USAGE:
%   [c,Ex]=itodensity(model,basis,cv);
% INPUTS
%   model  : a structure variable describing model (see below)
%   basis : a function definition structure specifying the
%             family of approximating functions
%   cv     : coefficient for the value function in control problems (optional)
% OUTPUTS
%   c      : approximation coefficients
%   Ex     : expected value associated with density 
%
% Note: no checks are performed to ensure that the process admits a
%       proper density.
%
% The model structure should have the following fields
%   func   : the name of the model function file
%   params : a cell array of model parameters
%
% If cv is passed, the model function file should follow the format 
%   of the model function files used by SCSOLVE.
% Otherwise the model function file should have the following syntax:
%   function out=funcfile(flag,s,additional parameters)
%   switch flag
%   case 'mu'
%     out = 
%   case 'sigma'
%     out = 
%   end
% The additional parameters must be in the same order as passed
% in the cell array model.params.

% Notes: Forms an approximation on a bounded interval using
%   p(x)=k*exp(int^x(2*mu(z)/sigma^2(z))dz)/sigma^2(x)
% where k is a constant that makes p(x) integrate to 1.

% Copyright (c) 1997-2002, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function [c,Ex]=itodensity(model,basis,cv)

 func=model.func;
 params=model.params;
 s=funnode(basis);

 if nargin<3                   % uncontrolled problem
   mu    = feval(func,'mu',s,params{:});
   sigma = feval(func,'sigma',s,params{:});
 else                          % controlled problem
   Vs    = funeval(cv,basis,s,1);
   x     = feval(func,'x',s,[],Vs,params{:});
   mu    = feval(func,'g',s,x,[],params{:});
   sigma = feval(func,'sigma',s,x,[],params{:});
 end
 sigma=sigma.*sigma;

 % Fit approximation to mu/sigma^2 and integrate
 c=funfitxy(basis,s,mu./sigma);
 temp=2*funeval(c,basis,s,-1);
 temp=temp-max(temp);                % normalize to avoid overflow

% Fit approximation to kernel (p)
 p=exp(temp)./sigma;
 c=funfitxy(basis,s,p);
 temp=funeval(c,basis,basis.b,-1); % determine constant
 c=c./temp;

% Compute expected value (if desired)
 if nargout>1
   p=p/temp;
   cc=funfitxy(basis,s,s.*p);
   Ex=funeval(cc,basis,basis.b,-1);
 end