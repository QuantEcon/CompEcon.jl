% LCPSTEP - Newton step for Array Linear Complementarity Problem

% Copyright (c) 1997-2010, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function [F,dx] = lcpstep(method,x,xl,xu,F,Fx)
if method(1)=='m'
  if nargout==1
    F = min(max(F,xl-x),xu-x);
  else
    F = min(max(F,xl-x),xu-x);
    dx = -arrayinvb(x,xl,xu,F,Fx);
    dx = min(max(dx,xl-x),xu-x);
  end
else
  if nargout==1
    F = arrayss(x,xl,xu,F);
  else
    [F,Fx] = arrayss(x,xl,xu,F,Fx);
    dx = -arrayinv(F,Fx);
    dx = min(max(dx,xl-x),xu-x);
  end
end

function y = arrayinvb(x,xl,xu,F,Fx)
[m,p,q] = size(Fx);
y = zeros(m,p);
AA = -eye(p);
for i=1:m
  A  = reshape(Fx(i,:,:),p,p);
  b  = F(i,:)';
  bl = (xl(i,:)-x(i,:))';
  bu = (xu(i,:)-x(i,:))';
  ind1 = b<=bl;
  ind2 = b>=bu;
  b(ind1)   = bl(ind1);
  A(ind1,:) = AA(ind1,:);
  b(ind2)   = bu(ind2);
  A(ind2,:) = AA(ind2,:);
  y(i,:) = (A\b)';
end

function y = arrayinv(F,Fx)
[m,p] = size(F);
y = zeros(m,p);
for i=1:m
  y(i,:) = (reshape(Fx(i,:,:),p,p)\reshape(F(i,:),p,1))';
end

function [Fnew,Fxnew] = arrayss(x,xl,xu,F,Fx)
if nargout<2
  Fnew = arrayssx(x,xl,xu,F);
else                            % Compute the Jacobian
  [Fnew,ff,aa] = arrayssx(x,xl,xu,F);
  [m,n] = size(x);
  Fxnew = repmat(ff,1,n).*reshape(Fx,m,n*n);
  % index for diagonal elements of the Jacobian
  ind = (0:n:n*(n-1)) + (1:n);
  Fxnew(:,ind) = Fxnew(:,ind)-aa;
end

function [Fnew,ffout,aaout] = arrayssx(x,xl,xu,F)
n = numel(x);
Fnew = zeros(size(x));
if nargout>1, ffout = zeros(size(x));
  if nargout>2, aaout = zeros(size(x)); end
end
if length(xl)==1, xl = xl+zeros(size(x)); end
if length(xu)==1, xu = xu+zeros(size(x)); end
for j = 1:n
  % compute phi+
  if isinf(xl(j)), d = F(j);
  else
    dxl = xl(j)-x(j);
    if abs(F(j))>abs(dxl), y = F(j); z = dxl;
    else                   y = dxl;  z = F(j);
    end
    z = z/y;
    dplus = sqrt(1+z*z);
    if y>0, d = y*(1+dplus+z);
    else    d = y*(z-((1-dplus)*(1-dplus)+z*z)/dplus/2);
    end
  end
  % compute phi-
  if isinf(xu(j)), Fnew(j) = d;
  else
    dxu = xu(j)-x(j);
    if abs(d)>abs(dxu), g = d;   h = dxu;
    else                g = dxu; h = d;
    end
    h = h/g;
    dminus = sqrt(1+h*h);
    if (g<0), Fnew(j) = g*(1+dminus+h);
    else      Fnew(j) = g*(h-((1-dminus)*(1-dminus)+h*h)/dminus/2);
    end
  end
  % compute Jacobian factors if requested
  if nargout>1
    if isinf(xu(j))
      ff = 1; aa = 1; bb = 0;
    else
      if g<0, dminus = -dminus; end
      temp1 = 1-1/dminus;  temp2 = 1-h/dminus;
      if abs(d)>abs(dxu), ff = temp1; aa = temp1; bb = temp2;
      else                ff = temp2; aa = temp2; bb = temp1;
      end
    end
    if isinf(xl(j)), aa = 0;
    else
      if y<0, dplus = -dplus; end
      temp1 = 1+1/dplus; temp2 = 1+z/dplus;
      if abs(F(j))>abs(dxl), ff = ff*temp1; aa = aa*temp2;
      else                   ff = ff*temp2; aa = aa*temp1;
      end
    end
    ffout(j) = ff;
    if nargout>2, aaout(j) = aa+bb; end
  end
end