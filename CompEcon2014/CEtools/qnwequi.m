%% QNWEQUI
%
%  Generates equidistributed nodes for computing the definite integral of
%  a real-valued function defined on a hypercube [a,b] in R^d.
%
%  Usage
%    [x,w] = qnwequi(n,a,b,type)
%  Input
%    n   : total number of nodes
%    a   : 1.d left endpoints
%    b   : 1.d right endpoints
%    type: type of sequence
%           N - Neiderreiter (default)
%           W - Weyl
%           H - Haber
%           R - pseudo Random
%  Output
%    x   : n.d quadrature nodes
%    w   : n.1 quadrature weights
%  Note
%    To compute definte integral of a real-valued function f defined on a
%    hypercube [a,b] in R^d, write a Matlab function f that returns an m.1
%    vector when passed an m.n matrix, and write [x,w]=qnwequi(n,a,b,*);
%    Intf=w'*f(x).

%  Copyright(c) 1997-2010
%   Mario J. Miranda - miranda.4@osu.edu
%   Paul L. Fackler  - paul_fackler@ncsu.edu

function [x,w] = qnwequi(n,a,b,type)

if any(a>b)
  error('In qnwequi: right endpoints must exceed left endpoints.')
end

global equidist_pp

if isempty(equidist_pp)
  equidist_pp = sqrt(primes(7920));   % good for d<=1000
end

d = max(length(n),max(length(a),length(b)));
n = prod(n);
if nargin<4
  type = 'N';
end

i = (1:n)';
switch upper(type(1))
  case 'N'                 % Neiderreiter
    j = 2.^((1:d)/(d+1));
    x = i*j;
    x = x-fix(x);
  case 'W'                 % Weyl
    j = equidist_pp(1:d);
    x = i*j;
    x = x-fix(x);
  case 'H'                 % Haber
    j = equidist_pp(1:d);
    x = (i.*(i+1)./2)*j;
    x = x-fix(x);
  case 'R'                 % pseudo-random
    x = rand(n,d);
  otherwise
    error('In qnwequi: unknown sequence requested.')
end

u = ones(n,1);
r = b-a;
x = a(u,:) + x.*r(u,:);
w = (prod(r)/n)*u;
