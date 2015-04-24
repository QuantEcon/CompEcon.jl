% FUNBCONV Converts among basis structure formats
% A function for converting basis structures to other formats.
% Conversion can only go "up", i.e., 'tensor' format can be
% converted to 'direct' or 'expanded' and 'direct' can be converted to
% 'expanded'
%
% The default is to convert to 'expanded'.
% The default order is the function itself (order 0).
%
% USES: chkfields, dprod
%
% See also: FUNBASEX

% Copyright (c) 1997-2000, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function B=funbconv(b,order,format)

if chkfields(b,{'vals';'order';'format'})
  disp('Basis structure has an invalid set of fields')
  error(' ')
end

d=size(b.order,2);
if nargin<3 || isempty(format), format='expanded'; end
if nargin<2 || isempty(order);
  if strcmp(format,'expanded')
    order=zeros(1,d);
  else
    order=b.order;
  end
end;

[numbas,d1]=size(order);
if d1~=d
  error('ORDER incompatible with basis structure');
end

if any(min(order,[],1)<b.order)% | any(max(order)>b.maxorder)
  error('Order of derivative operators exceeds those of basis');
end

% do this if you want to end up in expanded form (from tensor or direct)
if strcmp(format,'expanded')
  B.vals=cell(numbas,1);
  B.order=order;
  B.format=format;

  if strcmp(b.format,'tensor')
    m=1;n=1;
    for j=1:d
      m=m*size(b.vals{1,j},1);
      n=n*size(b.vals{1,j},2);
    end
    for i=1:numbas
      B.vals{i}=b.vals{order(i,d)-b.order(d)+1,d};
      for j=d-1:-1:1
        B.vals{i}=kron(B.vals{i},b.vals{order(i,j)-b.order(j)+1,j});
      end
    end
  elseif strcmp(b.format,'direct')
    n=1;for j=1:d, n=n*size(b.vals{1,j},2); end
    for i=1:numbas
      B.vals{i}=b.vals{order(i,d)-b.order(d)+1,d};
      for j=d-1:-1:1
        B.vals{i}=dprod(B.vals{i},b.vals{order(i,j)-b.order(j)+1,j});
      end
    end
  elseif strcmp(b.format,'expanded')
    %disp('WARNING: Basis is already in expanded form')
    B=b;
  else
    error('Improper basis format')
  end

% do this if you want to end up in direct form (from tensor)
elseif strcmp(format,'direct') && strcmp(b.format,'tensor')
  B.vals=cell(numbas,d);
  B.order=order;
  B.format=format;
  ind=cell(1,d);
  for j=1:d
    for i=1:size(b.vals,1)
      if ~isempty(b.vals{i,j})
        ind{j}=(1:size(b.vals{i,j},1))';
        break
      end
    end
  end
  [ind{:}]=gridmake(ind);
  for j=1:d
    for i=1:numbas
      if ~isempty(b.vals{i,j})
        B.vals{i,j}=b.vals{i,j}(ind{j},:);
      end
    end
  end
else
  disp('Not implemented for this option')
  error(' ')
end


% CHKFIELDS Checks if a variable S is a valid structure with fields F
% USAGE
%   errcode=chkfields(s,f);
% Returns an error code with values:
%   0 : no errors
%   1 : s is not a structure
%   2 : The number of fields in s differs from the number of elements of F
%   3 : The field list of s does not match that of F
%
% NOTE: to check fields interactively use fields(s)

% Copyright (c) 1997-2000, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function errcode=chkfields(s,f);

errcode=0;
if ~isstruct(s)
  errcode=1;
else
  ff=fieldnames(s);
  if any(size(ff)~=size(f)) errcode=2;
  elseif ~all(strcmp(ff,f)) errcode=3;
  end
end
