% ICRES Residual function for impulse control solver
% Used by ICSOLVE

% Copyright (c) 1997-2002, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function [e,c]=icres(y,x,func,params,xindex,basis,z, ...
                      Phi0,Phi1,Phi2,phi0,phi1,phi2,F)
  x(xindex>0) = y; 
  n = size(z,1); 
  a = x(1,1);
  r = x(2,1)-a;
  S = r*z+a;

  rho   = feval(func,'rho',S,params{:});
  mu    = feval(func,'mu',S,params{:});
  H     = spdiags(rho,0,n,n)*Phi0-spdiags(mu/r,0,n,n)*Phi1;
  sigma = feval(func,'sigma',S,params{:});
  if ~isempty(sigma)
     H = H - spdiags(sigma.*sigma/(2*r^2),0,n,n)*Phi2;
  end
  f = feval(func,'f',S,params{:});
  
  G = []; g = [];
  if xindex(1,1)>=0
    if F(1)==0  
      R = feval(func,'R+',[x(1) x(1)],params{:});
      H = [H;phi1(1,:)];
      f = [f;r*R(2)];
      if xindex(1,1)>0, G = [G;phi2(1,:)]; g = [g;r.^2*R(4)]; end
    else
      R = feval(func,'R+',x(1,:),params{:});
      H = [H;phi0(1,:)-funbase(basis,(x(1,2)-x(1,1))/r)];
      f = [f;R(1)-F(1)];
      if xindex(1,1)>0, G = [G;phi1(1,:)]; g = [g;r*R(2)]; end
      if xindex(1,2)>0, G = [G;funbase(basis,(x(1,2)-x(1,1))/r,1)]; g = [g;-r*R(3)]; end
    end
  end
  if xindex(2,1)>=0
     if F(2)==0
      R = feval(func,'R-',[x(2) x(2)],params{:});
      H = [H;phi1(2,:)];
      f = [f;r*R(2)];
      if xindex(2,1)>0, G = [G;phi2(2,:)]; g = [g;r.^2*R(4)]; end
    else
      R = feval(func,'R-',x(2,:),params{:});
      H = [H;phi0(2,:)-funbase(basis,(x(2,2)-x(1,1))/r)];
      f = [f;R(1)-F(2)];
      if xindex(2,1)>0, G = [G;phi1(2,:)]; g = [g;r*R(2)]; end
      if xindex(2,2)>0, G = [G;funbase(basis,(x(2,2)-x(1,1))/r,1)]; g = [g;-r*R(3)]; end
    end
  end

  c = H\f;
  e = G*c-g;
  