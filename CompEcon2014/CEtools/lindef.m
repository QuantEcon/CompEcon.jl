% LINDEF Computes standard breakpoints for linear spline
% USAGE
%  [n,a,b,parms]=LINDEF(breaks,evennum);
%  [n,a,b,parms]=LINDEF([a,b],num);
% INPUTS
%   breaks  : a breakpoint sequence
%   evennum : 1 if breakpoints are evenly spaced, 0 otherwise
% or
%   [a,b]   : left and right endpoints of approximation interval
%   num     : number of evenly spaced breakpoints
%
% Standard breakpoints are evenly spaced
%
% See also: FUNDEF

% Copyright (c) 1997-2000, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function [n,a,b,parms]=lindef(breaks,evennum);

if nargin<1 | isempty(breaks)
  error('Breakpoint information must be passed')
end
if nargin<2
  evennum=0;
end

n=length(breaks);
if n<2
  error('breakpoint sequence must contain at least two elements')
end
if any(diff(breaks)<=0)
  error('breakpoint sequence must be increasing');
end
if evennum~=0 
  if length(breaks)==2
    breaks=linspace(breaks(1),breaks(2),evennum)'; 
  else 
    if length(breaks)<2
      error('breakpoint must have at least two elements')
    end
    if any(abs(diff(diff((breaks))))>5e-15*mean(abs(breaks)))
      error('breakpoint sequence is not evenly spaced')
    end
    evennum=length(breaks);
  end
end

n=length(breaks);
a=breaks(1);
b=breaks(end);
parms={breaks,evennum};


