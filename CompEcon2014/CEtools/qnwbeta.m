%% QNWBETA
%
%  Generates Gaussian quadrature nodes and probability weights for Beta 
%  distribution with shape parameters a and b.
%
%  Usage
%    [x,w] = qnwbeta(n,a,b)
%  Input
%    n   : number of nodes
%    a   : Shape parameter (must be positive, default=1)
%    b   : Shape parameter (must be positive, default=1)
%  Output
%    x   : n.1 quadrature nodes
%    w   : n.1 probability weights
%  Note
%    Beta distribution is defined on interval [0,1].  Its mean and variance 
%    are mu=a/(a+b) and var=ab/((a+b)^2(a+b+1)).  Beta(1,1) is the uniform 
%    distribution.  To replicate mean mu and standard deviation std, set
%    a = (1-mu)*mu^2/std^2-mu and b = a*(1/mu-1). 
%  Note
%    To compute Ef(X) when f is real-valued and X is Beta(a,b), write a 
%    Matlab function f that returns m.1 vector when passed an m.1 vector, 
%    and write [x,w]=qnwbeta(n,a,b); Ef=w'*f(x).
%  Note
%    Based on an algorithm in Press, Teukolsky, Vetterling, and Flannery,
%    "Numerical Recipes in FORTRAN", 2nd ed. Cambridge U. Press, 1992.

%  Copyright(c) 1997-2014
%   Mario J. Miranda - miranda.4@osu.edu
%   Paul L. Fackler  - paul_fackler@ncsu.edu

function [x,w] = qnwbeta(n,a,b)

if nargin<2, a=1; end
if nargin<3, b=1; end
if a<0, error('In qnwbeta: a must be positive.'), end
if b<0, error('In qnwbeta: b must be positive.'), end

a = a-1;
b = b-1;
maxit = 25;
x = zeros(n,1);
w = zeros(n,1);
for i=1:n
  % Reasonable starting values
  switch i
    case 1
      an = a/n;
      bn = b/n;
      r1 = (1+a)*(2.78/(4+n*n)+0.768*an/n);
      r2 = 1+1.48*an+0.96*bn+0.452*an*an+0.83*an*bn;
      z = 1-r1/r2;
    case 2
      r1 = (4.1+a)/((1+a)*(1+0.156*a));
      r2 = 1+0.06*(n-8)*(1+0.12*a)/n;
      r3 = 1+0.012*b*(1+0.25*abs(a))/n;
      z = z-(1-z)*r1*r2*r3;
    case 3
      r1 = (1.67+0.28*a)/(1+0.37*a);
      r2 = 1+0.22*(n-8)/n;
      r3 = 1+8*b/((6.28+b)*n*n);
      z = z-(x(1)-z)*r1*r2*r3;
    case n-1
      r1 = (1+0.235*b)/(0.766+0.119*b);
      r2 = 1/(1+0.639*(n-4)/(1+0.71*(n-4)));
      r3 = 1/(1+20*a/((7.5+a)*n*n));
      z = z+(z-x(n-3))*r1*r2*r3;
    case n
      r1 = (1+0.37*b)/(1.67+0.28*b);
      r2 = 1/(1+0.22*(n-8)/n);
      r3 = 1/(1+8*a/((6.28+a)*n*n));
      z = z+(z-x(n-2))*r1*r2*r3;
    otherwise
      z = 3*x(i-1)-3*x(i-2)+x(i-3);
  end
  ab = a+b;
  % Rootfinding iterations
  for its=1:maxit
    temp = 2+ab;
    p1 = (a-b+temp*z)/2;
    p2 = 1;
    for j=2:n
      p3 = p2;
      p2 = p1;
      temp = 2*j+ab;
      aa = 2*j*(j+ab)*(temp-2);
      bb = (temp-1)*(a*a-b*b+temp*(temp-2)*z);
      c = 2*(j-1+a)*(j-1+b)*temp;
      p1 = (bb*p2-c*p3)/aa;
    end
    pp = (n*(a-b-temp*z)*p1+2*(n+a)*(n+b)*p2)/(temp*(1-z*z));
    z1 = z;
    z = z1-p1./pp;
    if abs(z-z1)<3e-14
      break
    end
  end
  if its>=maxit
    error('In qnwbeta: failure to converge.')
  end
  x(i) = z;
  w(i) = temp/(pp*p2);
end
x = (1-x)/2;
w = w*exp(gammaln(a+n)+gammaln(b+n)-gammaln(n+1)-gammaln(n+ab+1));
w = w/(2*exp(gammaln(a+1)+gammaln(b+1)-gammaln(ab+2)));