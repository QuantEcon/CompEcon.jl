%% FUNDEFN 
%
%   Defines a basis function family.
%
% Usage
%   basis = fundefn(basistype,n,a,b,order)
% Let
%   d = dimension of the basis function domain
% Input
%   basistype : string indicating basis function type ('cheb','spli' or 'lin')
%   n         : 1.d vector of integers > 1 indicating order of approximation per dimension
%   a         : 1.d vector of left endpoints of approximation intervals per dimension
%   b         : 1.d vector of right endpoints of approximation intervals per dimension
%   order     : for 'spli' basistype, the order of the spline (default: 3 for cubic)
% Output
%   basis  : a function family structure
%
% Uses: fundef
% See : fundef, funeval, funbase

% Copyright(c) 1997-2010
%   Mario J. Miranda - miranda.4@osu.edu
%   Paul L. Fackler  - paul_fackler@ncsu.edu

function basis = fundefn(basistype,n,a,b,order)

d = length(n);
if length(a) ~= d, error('a must be same dimension as n'); end
if length(b) ~= d, error('a must be same dimension as n'); end
if any(a>b), error('left endpoint must be less than right endpoint'); end
if any(n<2), error('n(i) must be greater than one'); end
if nargin<5 || isempty(order), order=3; end

params = cell(1,d);
switch basistype
  case 'cheb', for i=1:d; params(i)= {{'cheb',n(i),a(i),b(i)}}; end
  case 'spli', for i=1:d; params(i)= {{'spli',[a(i);b(i)],n(i)-order+1,order}}; end
  case 'lin',  for i=1:d; params(i)= {{'lin',[a(i);b(i)],n(i)}}; end
  otherwise,   error ('basis type must be ''cheb'',''spli'', or ''lin''')
end
basis = fundef(params{:});