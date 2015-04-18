%% FDJAC
%
%  Computes two-sided finite difference Jacobian
%
%  Usage
%    J = fdjac(f,x,varargin)
%  Input
%    f        : name of function of form fval = f(x,varargin)
%    x        : evaluation point
%    varargin : parameters passed to function f
%  Output
%    J        : finite difference Jacobian
%  Options
%    These options may be set by user with OPTSET (defaults in parentheses):
%    tol      : factor used to set step size (eps^(1/3))

%  Copyright(c) 1997-2010
%   Mario J. Miranda - miranda.4@osu.edu
%   Paul L. Fackler  - paul_fackler@ncsu.edu

function J = fdjac(f,x,varargin)

% Set options to defaults, if not set by user with OPTSET (see above)
tol = optget(mfilename,'tol',eps.^(1/3));

% Compute stepsize
h = tol.*max(abs(x),1);
xh1 = x+h; 
xh0 = x-h;
h   = xh1-xh0;

% Compute finite difference Jacobian
for j=1:length(x);
  xx = x;
  xx(j) = xh1(j); f1=feval(f,xx,varargin{:});
  xx(j) = xh0(j); f0=feval(f,xx,varargin{:});
  J(:,j) = (f1-f0)/h(j);
end

% One-sided derivative
% % Compute stepsize 
% h = tol.*max(abs(x),1);
% xh1 = x+h; 
% xh0 = x;
% h   = xh1-xh0;
% 
% % Compute finite difference Jacobian
% for j=1:length(x);
%   xx = x;
%   xx(j) = xh1(j); f1=feval(f,xx,varargin{:});
%   xx(j) = xh0(j); f0=feval(f,xx,varargin{:});
%   J(:,j) = (f1-f0)/h(j);
% end