% SPLIBASE Computes polynomial spline basis.
% USAGE
%   [B,x,k,breaks]=splibase(breaks,evennum,k,x,order);
% INPUTS
%   breaks  : user specified breakpoint sequence
%             (default: evenly spaced non-repeated breakpoints)
%   evennum : non-zero if breakpoints are all even
%   k       : polynomial order of the spline's pieces (default: 3, cubic)
%   x       : vector of the evaluation points (default: k-point averages of breakpoints)
%   order   : the order of differentiation (default: 0)
%             if a vector, SPLIBASE returns a cell array
%             otherwise it returns a matrix
% OUTPUTS
%   B : a kxn basis matrix
%   x : evaluation points (useful if defaults values are computed)
%
% Note: The number of basis functions is n=length(breaks)+k-1
%
% USES: splidop
%
% See also: SPLINODE, SPLIDEF, SPLIDOP, FUNBASE

% Copyright (c) 1997-2000, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function [B,x]=splibase(breaks,evennum,k,x,order)

if nargin<3, error('At least three parameters must be passed'); end
if nargin<4, x=[]; end
if nargin<5 | isempty(order), order=0; end

% GET DEFAULTS
if isempty(k), k=3; end

if isempty(x)
  x=splinode(breaks,evennum,k);
end

% A FEW CHECKS
if k<0
  error(['Incorrect value for spline order (k): ' num2str(k)]);
end
if min(size(breaks))>1
  error('''breaks'' must be a vector');
end
if any(order>=k)
  error('Order of differentiation must be less than k');
end
if size(x,2)>1
  error('x must be a column vector')
end

p=length(breaks);
m=size(x,1);
minorder=min(order);

% Augment the breakpoint sequence
n=length(breaks)+k-1; a=breaks(1); b=breaks(end);
augbreaks=[a(ones(k-minorder,1));breaks(:);b(ones(k-minorder,1))];

% The following lines determine the maximum index of
%   the breakpoints that are less than or equal to x,
%   (if x=b use the index of the next to last breakpoint).
%  [temp,ind]=sort([-inf;breaks(2:end-1);x(:)]);
%  temp=find(ind>=p);
%  j=ind(temp)-(p-1);
%  ind=temp-(1:m)';
%  ind(j)=ind(:)+(k-minorder);    % add k-minorder for augmented sequence
ind=lookup(augbreaks,x,3);

% Recursively determine the values of a k-order basis matrix.
% This is placed in an (m x k+1-order) matrix
bas=zeros(m,k-minorder+1);
bas(:,1)=ones(m,1);
B=cell(length(order),1);
if max(order)>0, D=splidop(breaks,evennum,k,max(order)); end % Derivative op
if minorder<0, I=splidop(breaks,evennum,k,minorder); end     % Integral op
for j=1:k-minorder
  for jj=j:-1:1
    b0=augbreaks(ind+jj-j);
    b1=augbreaks(ind+jj);
    temp=bas(:,jj)./(b1-b0);
    bas(:,jj+1)=(x-b0).*temp+bas(:,jj+1);
    bas(:,jj)=(b1-x).*temp;
  end
  % as now contains the order j spline basis
  ii=find((k-j)==order);
  if ~isempty(ii)
    ii=ii(1);
    % Put values in appropriate columns of a sparse matrix
    r=(1:m)'; r=r(:,ones(k-order(ii)+1,1));
    c=(order(ii)-k:0)-(order(ii)-minorder);
    c=c(ones(m,1),:)+ind(:,ones(k-order(ii)+1,1));
    B{ii}=sparse(r,c,bas(:,1:k-order(ii)+1),m,n-order(ii));
    % If needed compute derivative or anti-derivative operator
    if order(ii)>0
      B{ii}=B{ii}*D{order(ii)};
    elseif order(ii)<0
      B{ii}=B{ii}*I{-order(ii)};
    end
    %B{ii}=full(B{ii});
  end
end

if length(order)==1, B=B{1}; end
