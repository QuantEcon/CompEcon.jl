%% FUNEVAL 
%
%  Computes values and derivatives of function y=f(x) defined as linear
%  combination of specified basis functions.
%
%  Usage
%    y = funeval(c,basis,x,order)
%  Let
%    dx = dimension of x (domain)
%    dy = dimension of y (range)
%    nx = number of x values input
%    n  = 1.dx number of basis functions per domain dimension
%    nc = number of basis functions and coefficients (prod(n))
%    p  = number of derivatives requested
%  Input
%    c      : nc.dy basis function coefficients
%    x      : nx.dx domain values at which function is evaluated or
%             1.dx cell array of coordinates in each domain dimension
%    basis  : dx-dimensional function basis
%    order  : p.dx derivatives requested (optional, default=zeros(1,dx))
%  Output
%    y      : nx.dy function values or nx.dy.p derivative values
%
%  Note: Optional input order defines orders of the differential operator
%  to be applied. For exaple, order = [0 1 1], returns the cross partial
%  with respect to the 2nd and 3rd arguments. If multiple evaluations are
%  desired, use a row for each.  For example, order = [0 0 0;1 0 0;0 1 0]
%  returns the function value and the first two partial derivatives.
%  Negative values in the order matrix produce integrals normalized to
%  equal 0 at the left endpoints of the approximation intervals.  If order
%  is not supplied, just the function value is evaluated.
%
%  Alternate Usage (not for general use)
%    [y,B] = funeval(c,basis,B,order)
%  Input (Other inputs as above)
%    B      : basis structure predifined with funbasex
%  Output (Other outputs as above)
%    B      : basis structure used to evaluate function
%
% funeval makes use of subfunctions funeval1, funeval2, and funeval3
% correspoding to basis storage formats 'tensor', 'direct' and 'expanded'.
% See also: fundefn, funbase, funbasex

% Copyright (c) 1997-2014, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function [y,B] = funeval(c,basis,B,order)

if nargin<1 || (~isa(basis,'struct') && ~isa(B,'struct'))
  error('Either basis coefficients or basis structure must be passed'); 
end
if nargin<3; B = []; end
if nargin<4, order = 0; end

if isempty(c)
  error('Missing basis coefficients')
end

if ~isa(B,'struct')   % B is not a basis structure so construct one
  d = basis.d;
  if size(B,2)~=d
    error(['x must have d=' num2str(d) ' columns'])
  end
  if size(order,2)==1
    order = order*ones(1,d); 
  end
  B = funbasex(basis,B,order);
else
  if size(order,2)==1, order = order*ones(1,size(B.order,2)); end
end
  
if isa(B.vals,'cell') 
  switch B.format                  % determine the evaluator to call
    case 'tensor'  
      y = funeval1(c,B,order);
    case 'direct'
      y = funeval2(c,B,order);
    case 'expanded'
      y = funeval3(c,B,order);
    otherwise 
      error('Basis structure has an invalid ''format'' field');
  end
else
  y = B.vals*c;
end


% FUNEVAL1 Evaluates a function using 'tensor' format
function f = funeval1(c,B,order)
if nargin<3 || isempty(order)        % evaluate the function only
  kk = 1;
  d = size(B.order,2);
  order = zeros(1,d);
else                                % evaluate according to order
  [kk,d] = size(order);
end
% reverse the order of evaluation: B(d)xB(d-1)x...xB(1)
order = fliplr(order+ones(size(order,1),1)*(size(B.vals,1)*(0:d-1)-B.order+1));
nx = 1; for j=1:d, nx = nx*size(B.vals{1,j},1); end
f = zeros(nx,size(c,2),kk);          % preallocate output matrix
for i=1:kk
  f(:,:,i) = ckronx(B.vals,c,order(i,:));
end

% FUNEVAL2 Evaluates function using 'direct' format
function f = funeval2(c,B,order)
if nargin<3 || isempty(order)        % evaluate the function only
  kk = 1;
  d = size(B.order,2);
  order = zeros(1,d);
else
  [kk,d] = size(order);
end
% reverse the order of evaluation: B(d)xB(d-1)x...xB(1)
order = fliplr(order+ones(size(order,1),1)*(size(B.vals,1)*(0:d-1)-B.order+1));
f = zeros(size(B.vals{1},1),size(c,2),kk);   % preallocate
for i=1:kk
  f(:,:,i) = cdprodx(B.vals,c,order(i,:));
end


% FUNEVAL3 Evaluates a function using 'expanded' format
function f = funeval3(c,B,order)
if nargin<3 || isempty(order)   % default: evaluate the function only
  if isa(B.vals,'cell')
    kk = length(B.vals);
    order = 1:kk;
  else
    kk = 1;
    order = 1;
  end
else                           % determine which bases are to be used
  kk = size(order,1);
end
if isa(B.vals,'cell')
  nx = size(B.vals{1},1);
  f = zeros(nx,size(c,2),kk);
  for i=1:kk
    % Determine which element of B.vals is the desired basis, if any
    ii = find(ismember(B.order,order(i,:),'rows'));
    if isempty(ii)
      error('Requested basis matrix is not available');
    end
    if length(ii)>1,
      warning('Redundant request in FUNEVAL3')
      ii = ii(1);
    end
    f(:,:,i) = B.vals{ii}*c;
  end
else
  nx = size(B.vals,1);
  f = zeros(nx,size(c,2),kk);
  for i=1:kk
    f(:,:,i) = B.vals*c;
  end
end