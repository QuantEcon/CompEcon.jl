% CKRONXI The product of repeated inverse Kronecker products and a matrix.
% USAGE
%   z=ckronxi(b,c)       Solves inv(B1xB2x...xBd)*c
%   z=ckronxi(b,c,ind)   Selects cell elements according to IND
% where x denotes Kronecker (tensor) product.
% The Bi are passed as a cell array B. 
% B must be a vector cell array containing 2-D numerical arrays.
% The dimensions of B and C must agree:
%   if size(c,1) must equal product over ind of size(b{ind},2).

% Copyright (c) 1997-2000, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function z=ckronxi(b,c,ind)

if nargin<2, error('At least two parameters must be passed'), end
if nargin<3, ind=1:length(b); end

if ~iscell(b)                         % b is a matrix: return b*c
  if size(b,2)~=size(c,1)
    error('b and c are not conformable')
  end
  z=b\c;
else                                  % b is a cell array
  d=length(ind);
  n=zeros(d,1);
  for i=1:d n(i)=size(b{ind(i)},2); end
  if prod(n)~=size(c,1)
    error('b and c are not conformable')
  end
  z=c';
  mm=1;
  for i=1:d
    m=prod(size(z))/n(i);
    z=reshape(z,m,n(i));
    z=b{ind(i)}\z';
    mm=mm*size(z,1);
  end
  z=reshape(z,mm,size(c,2));
end
