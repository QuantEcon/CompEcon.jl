%% FUNBASE
%
%  Computes a basis matrix
%
%  Usage
%    B = funbase(basis,x,order)
%  Input
%    basis  : a structure defining a family of basis functions (created by fundef)
%    x      : m.d matrix or a 1xd cell array of columns vectors (created by funnode)
%    order  : 1xd vector indicating order of differentiation (default: zeros(1,d)))
%  Output
%    B      : m.prod(basis.n) basis matrix
%  Uses 
%    funbasex

%  Copyright(c) 1997-2010
%   Mario J. Miranda - miranda.4@osu.edu
%   Paul L. Fackler  - paul_fackler@ncsu.edu


function B = funbase(basis,x,order)

if nargin<1 || ~isstruct(basis)
  error('In funbase: "basis" must be a structure');
end
if nargin<3 || isempty(order)
  order = zeros(1,basis.d);
end
if nargin<2 || isempty(x)
  x = funnode(basis);
end

if size(order,1)>1,
  warning('In funbase: "order" should have only one row')
  order = order(1,:);
end

B = funbasex(basis,x,order,'expanded');
B = B.vals{1};