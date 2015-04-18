%% SHOW 
%
%  Displays a matrix with specified numbers of decimal places
%
%  Usage
%    show(x,d)
%  Input
%    x   : m.n matrix of numbers
%    d   : 1.n vector indicating # of decimal places for each column of x
%  Note
%    Input d is optiona; if omitted, 0 is used for columns of integers and 
%    4 is used otherwise.  Columns are separated by two blank spaces.

%  Copyright(c) 1997-2010
%   Mario J. Miranda - miranda.4@osu.edu
%   Paul L. Fackler  - paul_fackler@ncsu.edu

function show(x,d)

space='  ';
x=squeeze(x);
if length(size(x))>2
  error('Not defined for multidimensional arrays')
end

x=full(x);
if nargin<2 || isempty(d)
  d=4*(any(x~=fix(x)));
end
n=size(x,2);
if size(d,2)<n
  if length(d)==1
    d=d+zeros(1,n);
  else
    d=[d zeros(1,n-size(d,2))];
  end
end

fmt=[];
for i=1:n
  if size(x,1)>1
    m=max(x(x(:,i)<inf & ~isnan(x(:,i)),i));
    s=min(x(x(:,i)>-inf & ~isnan(x(:,i)),i))<0;
  else
    m=x(1,i);
    s=x(1,i)<0;
  end
  m=fix(log(m+eps)./log(10))+1;
  temp=real(m+s+d(i)+(d(i)>0));
  fmt=[fmt  space '%' num2str(temp) '.' num2str(d(i)) 'f'];
end
fmt=[fmt '\n'];
fprintf(fmt,x')