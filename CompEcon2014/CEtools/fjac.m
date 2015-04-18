%% FJAC
%  
%  Computes two-sided finite difference Jacobian using a more flexible 
%  input/output arrangement than fdjac.
%
%  Usage
%    J = fjac(f,ind,x1,x2,...)
%  Input
%    f         : name of function of form fval = f(x1,x2,...)
%    ind       : index of input/output arguments for Jacobian
%    x1,x2,... : input arguments for f
%  Output
%    J         : finite difference Jacobian
%  Example: 
%    If f has calling syntax [y1,y2]=f(x1,x2,x3), then
%       J=fjac(f,[3,2],x1,x2,x3)
%    computes approximation for dy2/dx3. Size of output wil be 
%    prod(size(y2)) by prod(size(y3)).
%  Options
%    Options may be set by user with OPTSET (defaults in parentheses):
%    tol    : factor used to set step size (eps^(1/3))

%  Copyright(c) 1997-2010
%   Mario J. Miranda - miranda.4@osu.edu
%   Paul L. Fackler  - paul_fackler@ncsu.edu

function J = fjac(f,ind,varargin)

% Set options to defaults, if not set by user with OPTSET (see above)
tol = optget(mfilename,'tol',eps.^(1/3));

% Determine indices of input and output variables
switch length(ind)
  case 2
    in = ind(1); out = ind(2);
  case 1
    in = ind; out = 1;
  case 0
    in = 1;  out = 1;
  otherwise
    error('Incorrectly specified index variable (ind)')
end

x = varargin{in}(:);
h   = tol.*max(abs(x),1);
xh1 = x+h; 
xh0 = x-h;
h   = xh1-xh0;
fval = cell(out,1);
for j=1:length(x)
  varargin{in}(:) = x;
  varargin{in}(j) = xh1(j); [fval{:}] = feval(f,varargin{:}); f1 = fval{out};
  varargin{in}(j) = xh0(j); [fval{:}] = feval(f,varargin{:}); f0 = fval{out};
  if j==1
    J = zeros(size(f1(:),1),size(x(:),1)); 
  end
  J(:,j) = (f1(:)-f0(:))/h(j);
end