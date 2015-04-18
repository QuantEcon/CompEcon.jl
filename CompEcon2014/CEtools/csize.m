% CSIZE Returns dimension information for cell arrays.
% USAGE
%   d=csize(b);
%   [d,n]=csize(b);
% If b is a cell array then:
%   d is the number of matrices in b
%   m is a dx2 matrix of row and column dimensions
%     or m and n are dx1 vectors
% If b is not a cell array
%   d=[]; this can be used to test if b is a cell array 
%         (isempty(d) is true is b is not a cell array)
%   m=size(b) or m=size(b,1) and n=size(b,2)

% Copyright (c) 1997-2000, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function [d,m,n]=csize(b)

if iscell(b)
  d=length(b);
  m=zeros(d,2);
  for i=1:d
    m(i,:)=size(b{i});
  end
else
  d=[];
  m=size(b);
end

if nargout==0
  disp([(1:d)' m])
elseif nargout==3
  n=m(:,2);
  m=m(:,1);
end

