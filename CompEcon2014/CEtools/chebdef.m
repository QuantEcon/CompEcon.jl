% CHEBDEF Defines parameters for Chebyshev polynomial functions
% USAGE
%   [n,a,b,parms]=chebdef(n,a,b);
% INPUTS
%   n    : the number of basis functions (1 plus the polynomial order)
%   a    : the left endpoint
%   b    : the right endpoint
%
% See also: CHEBBASE, FUNDEF

% Copyright (c) 1997-2000, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function [n,a,b,parms]=chebdef(n,a,b);

if nargin<3, error('3 parameters must be specified'), end
if n<=0 | fix(n)~=n,
  error('n must be a positive integer')
end
if (a>=b)
  error('Left end-point (a) must be less than right end-point (b)');
end

parms={n,a,b};