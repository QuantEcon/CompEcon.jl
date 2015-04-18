% CHEBDOP  Creates differential operator matrices for Chebyshev polynomials.
% USAGE
%   D=chebdop(n,a,b,order)          operator for ORDERth derivative
%   [D,n,a,b,parms]=CHEBDOP(n,a,b,order)  
%        Arguments 2-5 contain altered parameters 
% Computes the matrix operator that maps a vector of Chebyshev
% coefficients into the vector of coefficients of the derivative.
% For order 1, element i,j of this matrix is equal to
%    2(j-1) if i<j and i+j is odd
% for the interval [-1,1].  For other intervals divide by (b-a)/2.
% Higher order derivatives operators are formed by taking multiple products.
% Example:
%    g(x)=B^{n}(x)*c;
%    g'(x)=B^{n-1}(x)*chebdop(n,a,b)*c;
%    g"(x)=B^{n-2}(x)*chebdop(n,a,b,2)*c;
% Negative values of ORDER produce an integral operator normalized
%   to be zero at the left endpoint (a).
%
% D is returned as a cell array of size abs(order) composed of 
%   sparse matrices.
%
% Orders 1 and 2 are stored as a global cell array to avoid recomputation
% This speeds up repeated calculations involving the first two derivatives
%
% Most users will not need this function; 
%   use CHEBBASE, FUNEVAL or FUNDOP instead.
% 
% See also: CHEBBASE, FUNEVAL, FUNDOP.

% Copyright (c) 1997-2000, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function [D,n,a,b,parms]=chebdop(n,a,b,order)

global cheb_dop 

  if nargin<3 error('3 parameters must be specified'), end
  if nargin<4; order=1; end
  
  if order>0
    % Use previously stored values for order=1 and order=2
    % this speeds up calculations considerably with repeated calls
    if ~isempty(cheb_dop) & order<=length(cheb_dop) & size(cheb_dop{1},2)==n
      D=cell(order,1);
      for ii=1:order
        D{ii}=cheb_dop{ii}*(1/((b-a).^ii));
      end
    else
      D=cell(max(2,order),1);
      j = 1:n; j = j(ones(n,1),:); i = j';
      [r,c] = find(rem(i+j,2)==1 & j>i);
      d = sparse(r,c,(4/(b-a))*(j(1,c)-1),n-1,n);
      d(1,:) = d(1,:)/2;
      D{1}=d;
      for ii=2:max(2,order)
        D{ii}=d(1:n-ii,1:n-ii+1)*D{ii-1};
      end
      % store values for order=1 and order=2
      cheb_dop{1}=(b-a)*D{1};
      cheb_dop{2}=((b-a).^2)*D{2};
      if order==1, D(2)=[]; end
    end
  elseif order<0
    D=cell(abs(order),1);
    nn=n-order;
    i=(0.25*(b-a))./(1:nn);
    d=sparse([1:nn 1:nn-2],[1:nn 3:nn],[i -i(1:nn-2)],nn,nn);
    d(1,1)=2*d(1,1);
    d0=((-1).^(0:nn-1)).*sum(d);
    D{1}=[d0(1:n);d(1:n,1:n)];
    for ii=-2:-1:order;
      D{-ii}=[d0(1:n-ii-1);d(1:n-ii-1,1:n-ii-1)]*D{-ii-1};
    end
  else
    D=speye(n);
  end

if nargout>1
  n=n-order;
  parms={n,a,b};
end