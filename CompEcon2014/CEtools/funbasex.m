%% FUNBASEX
%
%  Creates basis structures for function evaluation (Use funbase to obtain
%  single basis matrices)
%
%  Usage
%    B = funbasex(basis,x,order,bformat)
%  Input
%    basis   : a structure defining a family of basis functions (default: created by fundef)
%    x       : mxd matrix or 1xd cell array of columns vectors (default: created by funnode)
%    order   : 1.d matrix (default: zeros(1,d)))
%    bformat : a string: 'tensor', 'direct' or 'expanded'
%  Output
%    B : a basis structure (defined below)
%    x : the computed evaluation points if none are passed
%  Defaults for bformat
%    'expanded' if d=1, otherwise
%    'tensor'   if x is a cell array
%    'direct'   if x is a matrix
%  B-Structure
%    vals     : cell array containing basis data (see exception below)
%    format   : 'tensor', 'direct', 'expanded'
%    order    : orders of differentiation ('expanded' format) or smallest
%               orders of differentiation and number of bases ('tensor' and 'direct' formats)
%  Note
%    Order Determines the # of basis matrices created:
%      for 'tensor' and 'direct' formats order should be:
%           1xd if only a single basis matrix is needed in each dimension
%           2xd specifying the minimum and maximum orders in each dimension
%           kxd listing of all desired basis matrices (only min and max of order used)
%     for 'expanded' format:
%           kxd listing of all desired basis matrices
%    If d=1 the format will be 'expanded' regardless of bformat
%  Uses
%    gridmake, funnode, funbconv
%  See
%    fundef, funnode, funveval, chebas, splibase, linbase, fourbas

%  Copyright(c) 1997-2010
%   Mario J. Miranda - miranda.4@osu.edu
%   Paul L. Fackler  - paul_fackler@ncsu.edu

function [B,x] = funbasex(basis,x,order,bformat)

if nargin<1 || ~isstruct(basis)
  error('A coefficient structure must be specified')
end
if nargin<2
  x = [];
end
if nargin<3 || isempty(order)
  order = 0;
end
if nargin<4
  bformat = [];
end

% Determine the problem dimension
d = length(basis.n);

% Expand ORDER if it has a single columns
if d>1 && size(order,2)==1
  order = order*ones(1,d);
end

% Initialize basis structure
m = size(order,1);
if m>1
  minorder = min(order)+zeros(1,d);
  numbases = (max(order)-minorder)+1;
else
  minorder = order+zeros(1,d);
  numbases = ones(1,d);
end

B = struct('vals',[],'order',minorder,'format',bformat);
B.vals = cell(max(numbases),d);

if isempty(x)
  x = funnode(basis);
end

if isempty(bformat)
  if isa(x,'cell')
    bformat = 'tensor';
  else
    bformat = 'direct';
  end
end

if d>1
  if ~isa(x,'cell') && strcmp(bformat,'tensor')
    error('Must pass a cell array to form a tensor format basis structure')
  end
  if isa(x,'cell') && strcmp(bformat,'direct')
    % it would be more efficient in this case to
    % use the tensor form to compute the bases and then
    % to use indexing to expand to the direct form
    %      x=gridmake(x); % convert to grid for direct form
  end
end

% notice we will be doing either tensor or direct here, _not_ expanded
if isa(x,'cell')
  B.format = 'tensor';
else
  B.format = 'direct';
end

% Compute basis matrices
switch B.format
  case 'tensor'
    for j=1:d
      if (m>1)
        orderj = unique(order(:,j));
      else
        orderj = order(1,j);
      end
      if length(orderj)==1
        B.vals{1,j} = feval([basis.bastype{j} 'base'],basis.parms{j}{:},x{j},orderj);
      else
        B.vals(orderj-minorder(j)+1,j) = feval([basis.bastype{j} 'base'],basis.parms{j}{:},x{j},orderj);
      end
    end
  case 'direct'
    for j=1:d
      if (m>1)
        orderj = unique(order(:,j));
      else
        orderj = order(1,j);
      end
      if length(orderj)==1
        % notice need to do this feval + splice string names together to get
        % function names. multiple dispatch takes care of this for us -- we just
        % need evalbase
        B.vals{1,j} = feval([basis.bastype{j} 'base'],basis.parms{j}{:},x(:,j),orderj);
      else
        B.vals(orderj-minorder(j)+1,j) = feval([basis.bastype{j} 'base'],basis.parms{j}{:},x(:,j),orderj);
      end
    end
end

% d=1, switch to expanded format
if size(B.vals,2)==1
  B.format = 'expanded';
  B.order = order;
  B.vals = B.vals(order+(1-min(order)));
  return
end

% Create expanded format
switch bformat
  case 'expanded'
    % funbconv converts from one basis representation to another. In julia we
     % have a function for that -- it's called `convert`
    B = funbconv(B,order,'expanded');
  case 'direct'
    if isa(x,'cell'), B = funbconv(B,order,'direct'); end
end
