%% ODESPLOT
%
%  Generates separatrix through saddle point of 2-dimensional ODE
%
%  ODESPX generates the separatrix of 2-dimensional first-order ODE
%    x'(t) = f(x(t))     t in [0,T]
%  by solving the ODE backwards in time, starting from near the saddle
%  point in the directions of the stable eigenvector. Here, x is 2.1
%  vector-valued function defined on time domain [0,T] and x' is its 2.1
%  vector-valued derivative with respect to time.
%
%  Usage
%   x = odespx(f,x,T,N,varargin)
%  Input
%    f        : name of velocity function (see below)
%    x        : presumed 2.1 saddle point
%    T        : time horizon in direction of each stable eigenvector
%    N        : number of time nodes in direction of each stable
%               eigenvector
%    varargin : parameters passed to velocity function (empty)
%  Output
%    x        : N.2 separatrix
%  NAN and infinite values removed from x before returning.

%  Copyright(c) 1997-2014
%    Mario J. Miranda - miranda.4@osu.edu

function odeplot(x,x1lim,x2lim)

N = size(x,1);
M = max(1,floor(N/100));
for j=M+1:N
  if mod(j,M)==0
    plot(x(j-M:j,1),x(j-M:j,2),'r')
    getframe;
    if x(j-M,1)<x1lim(1)||x(j-M,1)>x1lim(2)||x(j-M,2)<x2lim(1)||x(j-M,2)>x2lim(2), break, end
  end
end
plot(x(:,1),x(:,2),'r')