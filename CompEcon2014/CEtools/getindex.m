% GETINDEX Finds the index value of a point
% USAGE
%   i = getindex(s,S);
% INPUTS
%   s : a 1xd vector or pxd matrix
%   S : an nxd matrix
% OUTPUT
%   i : a p-vector of integers in {1,...,n} indicating the row of
%         S that most closely matches each row in s
%
% Example:
%  S=[0 0; 0 1; 1 0; 1 1]; s=[0 0; 0 0; 1 0; 0 0; 1 1];
%  getindex(s, S) returns
%  [1; 1; 3; 1; 4]
%
% Coded as a MEX file in C.

% Copyright (c) 1997-2002, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function i = getindex(s,S);
[n,d] = size(S);
p = size(s,1);
if p==1
  [junk,i] = min(sum(abs(s(ones(n,1),:)-S),2));
else
  [S,s]=gridmake(S,s);
  z=sum(abs(s-S),2);
  [junk,i]=min(reshape(z,n,p));
  i=i';
end