% FINDSTATE Calibrates an asset pricing model to data
% USAGE
%   s=findstate(c,basis,s,v,w);
% Given a model in which an asset price is 
% approximated by phi(s)c (funeval(c,basis,s))
% and observed prices are v, this function computes
% the weighted least squares estimate of s:
%   min (phi(s)c-v)diag(w)(phi(s)c-v)'
% where s is 1xd, phi(s) is 1xn, c is nxm and v and w are 1xm.
% w is optional; if omitted unweighted least squares is used.

% Copyright (c) 1997-2001, Mario J. Miranda & Paul L. Fackler
% miranda.4@osu.edu, paul_fackler@ncsu.edu

function s=findstate(c,basis,s,v,w)
maxit=100;
tol=sqrt(eps);

[n,m]=size(c);

if nargin>4 & ~isempty(w)
  c=c*spdiags(sqrt(w(:)),0,m,m);
  v=v*spdiags(sqrt(w(:)),0,m,m);
end

s=s(:)';
d=size(s,2);

for i=1:maxit
  s0=s;
  b=funbasex(basis,s,[0;1]);
  X=funeval(c,basis,b,eye(d));
  y=funeval(c,basis,b,0)-v;
  s=s-y/X;
  if all(abs(s-s0)<tol), break; end
end
