%% ODESPX
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

function x = odespx(f,x0,T,x1lim,x2lim,varargin)

N = 1000;

J = fdjac(f,x0,varargin{:});
[V,D] = eig(J);
i = find(diag(D)<0);
j = find(diag(D)>0);
if ~isreal(D)|i==0
  warning('x is not saddle point or stable steady-state')
  x = [];
  return
end

[~,i] = min(diag(D));
delx = 0.0001*V(:,i);
t = nodeunif(N,0,-T);
h = [0;diff(t)];

xsp = zeros(N,2);
x = x0+delx;
xsp(1,:) = x;
for i = 2:N
  hh = h(i);
  f1 = feval(f,x,varargin{:})*(hh/2);
  f2 = feval(f,x+f1,varargin{:})*hh;
  f3 = feval(f,x+f2/2,varargin{:})*hh;
  f4 = feval(f,x+f3,varargin{:})*(hh/2);
  x = x+(f1+f2+f3+f4)/3;
  xsp(i,:,:) = x;
  if x(1)<x1lim(1)||x(1)>x1lim(2)||x(2)<x2lim(1)||x(2)>x2lim(2), break, end
end
xsp(i+1:N,:) = [];
xsp = real(xsp);

xsn = zeros(N,2);
x = x0-delx;
xsn(1,:) = x;
for i = 2:N
  hh = h(i);
  f1 = feval(f,x,varargin{:})*(hh/2);
  f2 = feval(f,x+f1,varargin{:})*hh;
  f3 = feval(f,x+f2/2,varargin{:})*hh;
  f4 = feval(f,x+f3,varargin{:})*(hh/2);
  x = x+(f1+f2+f3+f4)/3;
  xsn(i,:,:) = x;
  if x(1)<x1lim(1)||x(1)>x1lim(2)||x(2)<x2lim(1)||x(2)>x2lim(2), break, end
end
xsn(i+1:N,:) = [];
xsn = real(xsn);

x = [xsn(end:-1:1,:);x0';xsp];
i = find(all(isfinite(x),2));
x = x(i,:);