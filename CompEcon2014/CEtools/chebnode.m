% CHEBNODE Computes standard nodes for Chebyshev polynomials
% USAGE
%   x=chebnode(n,a,b)
% Evaluates Gaussian nodes for the order n-1 Chebyshev polynomial
%   and transforms them to the interval [a,b].
%
% Setable options (use OPTSET):
%   nodetype   0 : Gaussian nodes (do not include endpoints)
%              1 : Gaussian nodes extended to endpoints
%              2 : Lobatto nodes  (includes endpoints)
%
% See also: CHEBBASE, FUNNODE.

% Copyright (c) 1997-2000, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function x=chebnode(n,a,b)

nodetype=optget('chebnode','nodetype',0);

s=(b-a)/2;
m=(b+a)/2;

if nodetype<2                           % usual nodes 
  k=pi*(0.5:(max(n)-0.5))';
  x=m(1)-cos(k(1:n(1))/n(1))*s(1);
  if nodetype==1                        % Extend nodes to endpoints
    aa=x(1);
    bb=x(end);
    x=(bb*a-aa*b)/(bb-aa)+(b-a)/(bb-aa)*x;
  end
else                                    % Lobatto nodes
  k=pi*(0:(max(n)-1))';
  x=m(1)-cos(k(1:n(1))/(n(1)-1))*s(1);
end