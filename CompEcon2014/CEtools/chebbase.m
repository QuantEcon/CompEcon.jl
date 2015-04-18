% CHEBBASE Computes basis matrices for Chebyshev polynomials
% USAGE
%   [B,x]=chebbase(n,a,b,x,order);
% INPUTS
%   n       : the number of basis functions (1 plus the polynomial order)
%   a       : the left endpoint
%   b       : the right endpoint
%   x       : k-vector of the evaluation points 
%             (default: roots of order n Chebyshev polynomial)
%   order   : the order of differentiation (default: 0)
%             if a vector, a cell array is returned
%             otherwise a matrix is returned
% OUTPUTS
%   B :  a kxn basis matrix or cell array of kxn basis matrices
%   x :  evaluation points (useful if defaults values are computed)
%
% USES: chebnode, chebdop
%
% See also: CHEBNODE, CHEBDOP, CHEBBASE, FUNBASE.

% Copyright (c) 1997-2001, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function [B,x]=chebbase(n,a,b,x,order)

  if nargin<3, error('3 parameters must be specified'), end
  if nargin<4, x=[]; end
  if nargin<5 || isempty(order), order=0; end

  minorder=min(0,min(order));

  if isempty(x); 
    x=chebnode(n,a,b);
    nodetype=optget('chebnode','nodetype');
  else
    nodetype=1;
  end
  
% Compute order 0 basis  
  if nodetype==0      % evaluate at standard nodes     
    temp=((n-0.5):-1:0.5)';
    bas=cos((pi/n)*temp*(0:(n-1-minorder)));
  else                % evaluate at arbitrary nodes
    bas=chebbasex(n-minorder,a,b,x);
  end
    
% Get bases for other orders  
  if length(order)==1
    if order~=0
      D=chebdop(n,a,b,order);
      B=bas(:,1:n-order)*D{abs(order)};
    else
      B=bas;
    end
  else
    B=cell(length(order),1);
    maxorder=max(order);
    if maxorder>0, D=chebdop(n,a,b,maxorder); end
    if minorder<0, I=chebdop(n,a,b,minorder); end
    for ii=1:length(order)
      if order(ii)==0
        B{ii}=bas(:,1:n);
      elseif order(ii)>0
        B{ii}=bas(:,1:n-order(ii))*D{order(ii)};
      else
        B{ii}=bas(:,1:n-order(ii))*I{-order(ii)};
      end
    end
  end