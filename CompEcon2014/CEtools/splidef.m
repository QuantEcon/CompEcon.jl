% SPLIDEF Defines default parameters for spline functions
% USAGE
%   [k,breaks]=splidef(breaks,evennum,k);
% INPUTS
%   breaks  : user specified breakpoint sequence
%             (default: evenly spaced non-repeated breakpoints)
%   evennum : non-zero if breakpoints are all even
%   k       : polynomial order of the spline's pieces (default: 3, cubic)
%
% See also:  SPLIBASE, FUNEVAL, FUNBASE

% Copyright (c) 1997-2000, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function [n,a,b,parms]=splidef(breaks,evennum,k);

if nargin<1
  error('At least one parameter must be passed');
end
if nargin<3 | isempty(k)
  k=3;                         % default is cubic splines
end
if nargin<2 | isempty(evennum)
  evennum=0;                         % default is unevenly spaced breakpoints
end

if prod(size(k,1))>1
  error(['Spline order (k) has improper size']);
end

if k<0
  error(['Spline order (k) is too small: ' num2str(k)]);
end

if length(breaks)<2
  error('breakpoint sequence must contain at least two elements');
end

if any(diff(breaks))<0
  error('Breakpoints must be non-decreasing');
end

if evennum==0
  if length(breaks)==2,  evennum=2;
  end
else
  if length(breaks)==2
    breaks=linspace(breaks(1),breaks(2),evennum)';
  else
    error('Breakpoint sequence must contain 2 values when evennum>0');
  end
end

n=length(breaks)+k-1;
a=breaks(1);
b=breaks(end);
parms={breaks,evennum,k};