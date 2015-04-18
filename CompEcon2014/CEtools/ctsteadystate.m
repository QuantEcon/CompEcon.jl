% CTSTEADYSTATE Finds deterministic steady state for continuous time models
% USAGE
%   [sstar,xstar,res]=ctsteadystate(model,basis,cv)
% INPUTS
%   model  : a structured variable (see scsolve for description)
%   basis : approximation function space structure variable (see fundef)
%   cv     : value function approximation coefficients
% OUTPUTS
%   sstar  : steady state state value
%   xstar  : steady state action value
%   res    : residual (can be used to check convergence)
%
% This function is a utility associated with scsolve for analyzing the
% behavior of continuous time stochastic control models

% Copyright (c) 1997-2002, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function [sstar,xstar,res]=ctsteadystate(model,basis,cv)
s=(basis.a+basis.b)/2;
sstar=broyden('ctssres',s,model,basis,cv);
if nargout>1
  xstar=feval(model.func,'x',sstar,[],funeval(cv,basis,sstar,1),model.params{:});
  if nargout>2
    res=ctssres(sstar,model,basis,cv);
  end
end