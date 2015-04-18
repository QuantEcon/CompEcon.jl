% NODEUNIF Computes uniform nodes for intervals R^n
% USAGE
%   [x,xcoord] = nodeunif(n,a,b);
%
% If dimension d of n, a and b is one, returns n by
% by 1 vector x containing n uniform nodes spanning
% the interval [a,b]. If dimension d>1, returns 1 by d
% cell array xcoord whose kth entry is the n(k) by 1
% vector of n(k) uniform nodes spanning the interval
% [a(k),b(k)]; also returns prod(n) by d matrix x of
% grid points created by forming Cartesian product of
% vectors in xcoord.

% Copyright (c) 1997-2000, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function [x,xcoord] = nodeunif(n,a,b)

d = length(n);
if d==1
   x = linspace(a,b,n)';
   xcoord = x;
else
   for k=1:d
      xcoord{k} = linspace(a(k),b(k),n(k))';
   end   
   x = gridmake(xcoord);
end
