% MINTERP Multidimensional linear interpolation
% USAGE
%   xinterp=minterp(scoord,x,s);
% INPUTS
%   scoord : a 1xd cell array of coordinate vectors
%              each must be sorted in ascending order
%            may be a vector rather than a cell array for d=1
%   x      : function values at the grid points defined
%              by scoord (use gridmake to expand gridpoints)
%              nxp array where n is the total number of grid
%              points defined by scoord (# of rows in gridmake(scoord))
%   s      : values at which to interpolate
%              s can be either:
%               1xd cell array of vectors
%                 to produce interpolates on a new grid
%               kxd matrix 
%                 to produce interpolates at k arbitrary points
%   evenspacing : 1 if all the vectors in scoord are evenly spaced
%                   (speeds up algorithm if true)
% OUTPUT
%   xinterp : interpolated values
%               Nxp if s is a cell array defining a grid with N points
%               kxp is s is a kxd matrix

% Copyright (c) 1997-2001, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function xinterp=minterp(scoord,x,s,evenspacing)

if nargin<4 |isempty(evenspacing)
  evenspacing=0;
end

d=size(s,2);
if d==1
  if isa(scoord,'cell')
    scoord=scoord{1};
  end
  if evenspacing, ind=evenlookup(scoord,s);
  else,           ind=lookup(scoord,s,3);
  end
  z=scoord(ind);
  xinterp=x(ind,:);
  ind=ind+1;
  z=(s-z)./(scoord(ind)-z);
  m=size(x,2);
  if m>1, z=z(:,ones(size(x,2),1)); end
  xinterp=xinterp.*(1-z)+x(ind,:).*z;
else
  if isa(s,'cell')       % use tensor products
    B=cell(1,d);
    for i=1:d
      B{i}=getbas(scoord{i},s{i},evenspacing);
    end
    xinterp=ckronx(B,x,d:-1:1);
  else                    % use direct (row-wise tensor) products
    B=cell(1,d);
    for i=1:d
      B{i}=getbas(scoord{i},s(:,i),evenspacing);
    end
    xinterp=cdprodx(B,x,d:-1:1);
  end
end

function B=getbas(scoord,s,evenspacing)
  m=size(s,1);
  n=size(scoord,1);
  if evenspacing, ind=evenlookup(scoord,s);
  else,           ind=lookup(scoord,s,3);
  end
  z=(s-scoord(ind))./(scoord(ind+1)-scoord(ind));
  B=sparse([1:m 1:m],[ind ind+1],[(1-z) z],m,n);
  

function ind=evenlookup(scoord,s)
  n=length(scoord);
  ind=ceil((s-scoord(1))*(n/(scoord(end)-scoord(1))));
  ind=min(max(ind,1),n-1);

