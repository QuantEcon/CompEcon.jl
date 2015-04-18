% SSMOOTH Reformulates an MCP as a semismooth function
% USAGE
%   [fxnew,Jnew]=ssmooth(x,a,b,fx,J);
% INPUTS
%    x : an evaluation point
%    a : lower bounds 
%    b : upper bounds 
%   fx : function values at x
%    J : Jacobian of f at x (optional, required if Jnew requested)
% OUTPUTS
%   fxnew : value of semi-smooth function at x
%    Jnew : Jacobian of function at x
%
% The reformulation uses
%   phi-(phi+(fx,a-x),b-x)
% where
%   phi+(y,z)=y+z+sqrt(y.^2+z^2)
%   phi-(y,z)=y+z-sqrt(y.^2+z^2)
%
% Uses: arrayssx

% Copyright (c) 1997-2000, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function [fxnew,Jnew]=ssmooth(x,a,b,fx,J)

try  % exists('arrayssx','file')
  if nargout<2
    fxnew=arrayssx(x,a,b,fx);
  else                            % Compute the Jacobian
    [fxnew,ff,xx]=arrayssx(x,a,b,fx);
    n=size(x,1);
    Jnew=spdiags(ff,0,n,n)*J-spdiags(xx,0,n,n);
  end
catch
  if length(a)==1,a=a+zeros(size(x)); end
  if length(b)==1,b=b+zeros(size(x)); end
  dainf = find(a==-inf);
  dbinf = find(b== inf);
  n=length(x);
  da = a-x;
  db = b-x;
  sq1=sqrt(fx.^2+da.^2);
  pval = fx+sq1+da;    pval(dainf)  = fx(dainf);
  sq2=sqrt(pval.^2+db.^2);
  fxnew = pval-sq2+db; fxnew(dbinf) = pval(dbinf);
  if nargout==2
    dpdy=1+fx./sq1;   dpdy(dainf)=1;  % y=fx, z=da
    dpdz=1+da./sq1;   dpdz(dainf)=0;
    dmdy=1-pval./sq2; dmdy(dbinf)=1;  % y=pval, z=db
    dmdz=1-db./sq2;   dmdz(dbinf)=0;
    ff=dmdy.*dpdy;                    % ff =  ds/df
    xx=dmdy.*dpdz+dmdz;               % xx = -ds/dx
    Jnew=spdiags(ff,0,n,n)*J-spdiags(xx,0,n,n);
  end
end