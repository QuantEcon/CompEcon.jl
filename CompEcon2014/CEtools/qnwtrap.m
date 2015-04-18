%% QNWTRAP
%
%  Generates Trapezoid rule quadrature nodes and weights for computing the
%  definite integral of a real-valued function defined on a hypercube [a,b]
%  in R^d.
%
%  Usage
%    [x,w] = qnwtrap(n,a,b)
%  Input
%    n   : 1.d number of nodes per dimension
%    a   : 1.d left endpoints
%    b   : 1.d right endpoints
%  Output
%    x   : prod(n).d quadrature nodes
%    w   : prod(n).1 quadrature weights
%  Note
%    To compute definte integral of a real-valued function f defined on a
%    hypercube [a,b] in R^d, write a Matlab function f that returns an m.1 
%    vector when passed an m.n matrix, and write [x,w]=qnwtrap(n,a,b); 
%    intf=w'*f(x).
%  Uses
%    ckron, gridmake, qnwtrap1

%  Copyright(c) 1997-2010
%   Mario J. Miranda - miranda.4@osu.edu
%   Paul L. Fackler  - paul_fackler@ncsu.edu

function [x,w] = qnwtrap(n,a,b)

d = length(n);
if nargin<2, a = zeros(1,d); end
if nargin<3, b = ones(1,d); end
if any(a>b)
  error('In qnwtrap: right endpoints must exceed left endpoints.') 
end

x = cell(1,d);
w = cell(1,d);
for i=1:d
   [x{i},w{i}] = qnwtrap1(n(i),a(i),b(i));
end
w = ckron(w(d:-1:1));
x = gridmake(x);


%% QNWTRAP1
%
%  Generates Trapezoid quadrature nodes and weights for computing the
%  definite integral of a real-valued function defined on an interval [a,b]
%  in R.
%
%  Usage
%    [x,w] = qnwtrap1(n,a,b)
%  Input
%    n   : number of nodes
%    a   : left endpoint
%    b   : right endpoint
%  Output
%    x   : n.1 quadrature nodes
%    w   : n.1 quadrature weights

function [x,w] = qnwtrap1(n,a,b)

if n<=1
  error('In qnwtrap: n must be integer greater than one.')
end

dx = (b-a)/(n-1);
x = (a:dx:b)';
w = dx*ones(n,1);
w([1;n]) = 0.5*w([1;n]);