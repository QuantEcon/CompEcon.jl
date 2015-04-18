% CKRON Repeated Kronecker products on a cell array of matrices.
% USAGE
%   z=ckron(b)     Solves (B1xB2x...xBd)
%   z=ckron(b,1)   Solves (inv(B1)xinv(B2)x...xinv(Bd))
% where x denotes Kronecker (tensor) product.
% The Bi are passed as a cell array B.

% Copyright (c) 1997-2000, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function z=ckron(b,invert)

if nargin<1, error('At least one parameter must be passed'), end
if nargin==1, invert=0; end

[d,m,n]=csize(b);
if invert & any(m~=n)
  error('Matrix elements must be square to invert');
end

if isempty(d)
  if invert z=inv(b); else z=b; end
else
  if invert z=inv(b{1})
  else z=b{1};
  end
  for i=2:d
    if invert
      z=kron(z,inv(b{i}));
    else
      z=kron(z,b{i});
    end
  end
end