%% QNWGAMMA
%
%  Generates Gaussian quadrature nodes and probability weights for Gamma
%  distribution with shape parameter a>0 and scale parameter b>0.  Also
%  generates Gaussian quadrature nodes and probability weights for
%  Chi-squared and exponential distributions, which are special cases of
%  the Gamma distribution (see below).
%
%  Usage
%    [x,w] = qnwgamma(n,a,b)
%  Input
%    n   : number of nodes
%    a   : shape parameter (must be positive, default=1)
%    b   : scale parameter (must be positive, default=1)
%  Output
%    x   : n.1 quadrature nodes
%    w   : n.1 probability weights
%  Note
%    Gamma distribution is defined on [0,inf).  Its mean and variance are
%    mu=ab and var=ab^2. Gamma(1,b) is the exponential distribution with
%    mean b.  Gamma(n/2,2) distribution is the Chi-squared distribution
%    with n degrees of freedom.
%  Note
%    To compute Ef(X) when f is real-valued and X is Gamma(a,b), write a
%    Matlab function f that returns m.1 vector when passed an m.1 vector,
%    and write [x,w]=qnwgamma(n,a,b); Ef=w'*f(x).
%  Note
%    Based on an algorithm in Press, Teukolsky, Vetterling, and Flannery,
%    "Numerical Recipes in FORTRAN", 2nd ed. Cambridge U. Press, 1992.

%  Copyright(c) 1997-2010
%   Mario J. Miranda - miranda.4@osu.edu
%   Paul L. Fackler  - paul_fackler@ncsu.edu

function [x,w] = qnwgamma(n,a,b)

if nargin<2, a=1; end
if nargin<3, b=1; end
if a<0, error('In qnwgamma: shape parameter a must be positive.'), end
if b<0, error('In qnwgamma: scale parameter b must be positive.'), end

a = a-1;
maxit = 10;
factor = -exp(gammaln(a+n)-gammaln(n)-gammaln(a+1));
x = zeros(n,1);
w = zeros(n,1);
for i=1:n
  % Reasonable starting values
  if i==1
    z = (1+a)*(3+0.92*a)/(1+2.4*n+1.8*a);
  elseif i==2
    z = z+(15+6.25*a)./(1+0.9*a+2.5*n);
  else
    j = i-2;
    z = z+((1+2.55*j)./(1.9*j)+1.26*j*a./(1+3.5*j))*(z-x(j))./(1+0.3*a);
  end
  % Rootfinding iterations
  for it=1:maxit
    p1 = 1;
    p2 = 0;
    for j=1:n
      p3 = p2;
      p2 = p1;
      p1 = ((2*j-1+a-z)*p2-(j-1+a)*p3)./j;
    end
    pp = (n*p1-(n+a)*p2)./z;
    z1 = z;
    z = z1-p1./pp;
    if abs(z-z1)<3e-14
      break;
    end
  end
  if it>=maxit
    error('In qnwgamma: failure to converge.')
  end
  x(i) = z;
  w(i) = factor/(pp*n*p2);
end
x = b*x;
