% FUNNODE
% Generates default interpolation nodes for a family of basis functions
%
% USAGE
%   [x,xcoord] = funnode(basis)
% INPUT
%   basis = a family of basis functions
% OUPUT
%   x      = prod(n) by d interpolation node grid
%   xcoord = if d>1, 1 by d cell array of node coordinates by dimension
% See also: FUNDEF, FUNBASE, FUNFITF, FUNFITXY, FUNEVAL.

% Copyright (c) 1997-2000, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function [x,xcoord] = funnode(basis)

d = basis.d;
if d==1 
  xcoord = feval([basis.bastype{1} 'node'],basis.parms{1}{:}); 
  x = xcoord;
else
  xcoord = cell(1,d);                        
  for j=1:d
    xcoord{j} = feval([basis.bastype{j} 'node'],basis.parms{j}{:});  
  end
  x = gridmake(xcoord);
end