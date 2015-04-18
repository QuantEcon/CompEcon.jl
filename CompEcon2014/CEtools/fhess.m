%% FHESS
%
%  Computes finite difference Hessian using a more flexible input/output 
%  arrangement than fdhess.
%
%  Usage
%    H = fhess(f,ind,x1,x2,...)
%  Input
%    f         : name of function of form fval = f(x1,x2,...)
%    ind       : index of input/output arguments for Hessian
%    x1,x2,... : input arguments for f
%  Output
%    H         : finite difference Hessian
%  Example: 
%    If f has calling syntax [y1,y2]=f(x1,x2,x3), then
%       H=fhess(f,[3,2],x1,x2,x3)
%    computes an approximation for d^2y2/dx3^2. Size of output will be
%    prod(size(y2)) by prod(size(x3)) by prod(size(x3)), except singleton
%    dimension will be elimiated if prod(size(y2))=1.
%  Options
%    Options may be set by user with OPTSET (defaults in parentheses):
%    tol      : factor used to set step size (eps^(1/4))
%    diagonly : computes just the diagonal elements of the Hessian (0)

%  Copyright(c) 1997-2010
%   Mario J. Miranda - miranda.4@osu.edu
%   Paul L. Fackler  - paul_fackler@ncsu.edu

function H = fhess(f,ind,varargin)

% Set options to defaults, if not set by user with OPTSET (see above)
tol      = optget(mfilename,'tol',eps.^(1/4));
diagonly = optget(mfilename,'diagonly',0);

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

fval = cell(out,1);
[fval{:}] = feval(f,varargin{:});
fx = fval{out}(:);
m = size(fx,1);

x = varargin{in}(:);
n = size(x,1);

% Compute stepsize
h  = tol.*max(abs(x),1);
xh = x+h;
h  = xh-x;
ee = sparse(1:n,1:n,h,n,n);

hh = h*h';
H  = zeros(m,n,n);

% Compute forward and backward steps
f1 = zeros(m,n);
f0 = zeros(m,n);
for i=1:n
  varargin{in}(:) = x+ee(:,i); [fval{:}] = feval(f,varargin{:}); f1(:,i) = fval{out}(:);
  varargin{in}(:) = x-ee(:,i); [fval{:}] = feval(f,varargin{:}); f0(:,i) = fval{out}(:);
  H(:,i,i) = (f1(:,i)+f0(:,i)-2*fx)./hh(i,i);
end

% Compute double steps
if ~diagonly
  for i=1:n
    for j=1:n
      if i~=j
        varargin{in}(:) = x+ee(:,i)-ee(:,j);
        [fval{:}] = feval(f,varargin{:}); fxx = fval{out}(:);
        H(:,i,j) = (f1(:,i)+f0(:,j)-fx-fxx)./hh(i,j);
      end
    end
  end
  H = (H+permute(H,[1 3 2]))./2; % transpose the 2 & 3 dimensions
end

if m==1 % if output is 1-D, return an nxn matrix
  H = squeeze(H);
end