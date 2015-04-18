% LINBASE Piecewise linear basis functions
% USAGE
%   [B,x]=linbase(breaks,evennum,x,order);
% INPUTS
%   breaks  : an nx1 vector of breakpoints
%   evennum : =n if breakpoints are all even, otherwise=0
%   x       : k-vector of the evaluation points 
%              (default: breaks)
%   order   : the order of differentiation (default: 0)
%             if a vector, SPLIBASE returns a cell array 
%             otherwise it returns a matrix
% OUTPUTS
%   B :  a kxn basis matrix or cell array of basis matrices
%   x :  evaluation points (useful if defaults values are computed) 
%
% See also: LINNODE, LINDEF, LINDOP, FUNBASE

% Copyright (c) 1997-2000, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function [B,x]=linbase(breaks,evennum,x,order)

if nargin<4 | isempty(order), order=0;
  if nargin<3, x=[];
    if nargin<2, evennum=0;
      if nargin<1, error('At least one parameter must be passed');
      end; end; end; end;

% GET DEFAULTS
if isempty(x)
  x=linnode(breaks,evennum);
end
n=length(breaks);

% If multiple orders are requested make recursive call
% Inefficient but easy to code!
k=length(order);
if k>1
  B=cell(k,1);
  for ii=1:k
    B{ii}=linbase(breaks,evennum,x,order(ii));
  end
  return
end

if order~=0               % recursively generate differential operators
  [D,n,a,b,parms]=lindop(breaks,evennum,order);
  B=linbase(parms{:},x)*D{end};
  return
end

m=size(x,1);

% Determine the maximum index of
%   the breakpoints that are less than or equal to x,
%   (if x=b use the index of the next to last breakpoint).
if evennum
  ind=fix((x-breaks(1)).*((n-1)./(breaks(end)-breaks(1))))+1;
  ind=min(max(ind,1),n-1);
else
  ind=lookup(breaks,x,3);
end

z=(x-breaks(ind))./(breaks(ind+1)-breaks(ind));
B=sparse([1:m 1:m],[ind ind+1],[(1-z) z],m,n);