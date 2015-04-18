%% QNWLEGE 
%
%  Generates Guass-Legendre quadrature nodes and weights for computing the
%  definite integral of a real-valued function defined on a hypercube [a,b]
%  in R^d.
%
%  Usage
%    [x,w] = qnwlege(n,a,b)
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
%    vector when passed an m.n matrix, and write [x,w]=qnwlege(n,a,b); 
%    intf=w'*f(x).
%  Uses
%    ckron, gridmake, qnwlege1

%  Copyright(c) 1997-2010
%   Mario J. Miranda - miranda.4@osu.edu
%   Paul L. Fackler  - paul_fackler@ncsu.edu

function [x,w] = qnwlege(n,a,b)

d = length(n);
if nargin<2, a = zeros(1,d); end
if nargin<3, b = ones(1,d); end
if any(a>b)
  error('In qnwlege: right endpoints must exceed left endpoints.') 
end

x = cell(1,d);
w = cell(1,d);
for i=1:d
   [x{i},w{i}] = qnwlege1(n(i),a(i),b(i));
end
w = ckron(w(d:-1:1));  % use reverse ordered tensor product
x = gridmake(x);


%% QNWLEGE1 
%
%  Generates Guass-Legendre quadrature nodes and weights for computing the
%  definite integral of a real-valued function defined on an interval [a,b]
%  in R.
%
%  Usage
%    [x,w] = qnwlege1(n,a,b)
%  Input
%    n   : number of nodes 
%    a   : left endpoints
%    b   : right endpoint
%  Output
%    x   : n.1 quadrature nodes
%    w   : n.1 quadrature weights
%  Note
%    Based on an algorithm in Press, Teukolsky, Vetterling, and Flannery,
%    "Numerical Recipes in FORTRAN", 2nd ed. Cambridge U. Press, 1992.

function [x,w] = qnwlege1(n,a,b)

if n<=1
  error('In qnwlege: n must be integer greater than one.')
end

maxit = 100;
m = fix((n+1)/2);
xm = 0.5*(b+a);
xl = 0.5*(b-a);
x = zeros(n,1);
w = x;
i = (1:m)';
z = cos(pi*(i-0.25)./(n+0.5));
for its=1:maxit
   p1 = 1;
   p2 = 0;
   for j=1:n
      p3 = p2;
      p2 = p1;
      p1 = ((2*j-1)*z.*p2-(j-1)*p3)./j;
   end
   pp = n*(z.*p1-p2)./(z.*z-1);
   z1 = z;
   z = z1-p1./pp;
   if abs(z-z1)<1e-14 
     break; 
   end
end
if its==maxit
   error('In qnwlege: failure to converge.')
end
x(i) = xm-xl*z;
x(n+1-i) = xm+xl*z;
w(i) = 2*xl./((1-z.*z).*pp.*pp);
w(n+1-i) = w(i);