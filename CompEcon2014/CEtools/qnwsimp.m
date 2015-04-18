%% QNWSIMP 
%
%  Generates Simpson's rule quadrature nodes and weights for computing the
%  definite integral of a real-valued function defined on a hypercube [a,b]
%  in R^d.
%
%  Usage
%    [x,w] = qnwsimp(n,a,b)
%  Input
%    n   : 1.d number of nodes per dimension (must be odd positive integers)
%    a   : 1.d left endpoints
%    b   : 1.d right endpoints
%  Output
%    x   : prod(n).d quadrature nodes
%    w   : prod(n).1 quadrature weights
%  Note
%    To compute definte integral of a real-valued function f defined on a
%    hypercube [a,b] in R^d, write a Matlab function f that returns an m.1 
%    vector when passed an m.n matrix, and write [x,w]=qnwsimp(n,a,b); 
%    intf=w'*f(x).
%  Uses
%    ckron, gridmake, qnwsimp1

%  Copyright(c) 1997-2010
%   Mario J. Miranda - miranda.4@osu.edu
%   Paul L. Fackler  - paul_fackler@ncsu.edu

function [x,w] = qnwsimp(n,a,b)

d = length(n);
if nargin<2, a = zeros(1,d); end
if nargin<3, b = ones(1,d); end
if any(a>b)
  error('In qnwsimp: right endpoints must exceed left endpoints.') 
end

x = cell(1,d);
w = cell(1,d);
for i=1:d
   [x{i},w{i}] = qnwsimp1(n(i),a(i),b(i));
end
w = ckron(w(d:-1:1));
x = gridmake(x);


%% QNWSIMP1
%
%  Generates Simpson's rule quadrature nodes and weights for computing the
%  definite integral of a real-valued function defined on an interval [a,b]
%  in R.
%
%  Usage
%    [x,w] = qnwsimp1(n,a,b)
%  Input
%    n   : number of nodes (must be an odd positive integer)
%    a   : left endpoint
%    b   : right endpoint
%  Output
%    x   : n.1 quadrature nodes
%    w   : n.1 quadrature weights

function [x,w] = qnwsimp1(n,a,b)

if n<=1
  error('In qnwsimp: n must be integer greater than one.')
end
if rem(n,2)==0
  warning('In qnwsimp: n must be odd integer - increasing by 1.')
  n = n+1;
end

dx = (b-a)/(n-1); 
x = (a:dx:b)';
w = reshape([2;4]*ones(1,(n+1)/2),n+1,1);
w = w(1:n); 
w(1) = 1; 
w(n) = 1;
w = (dx/3)*w;