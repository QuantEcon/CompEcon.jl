% LINDOP Differential operator for a piecewise linear function
% USAGE
%   D=lindop(breaks,evennum,order);
% INPUTS
%   breaks  : an nx1 vector of breakpoints
%   evennum : =n if breakpoints are all even, otherwise=0
%   order   : the desired order of differentiation (a scalar)
% OUTPUTS
%   D : a cell array of matrices of size abs(order) by 1.
%   n, a, b, parms : describes the characteristics of the altered
%                    family of functions (how fundef is altered)
%                    this is used by FUNDOP
%
% A piecewise linear function with n-1 pieces can be described
% by the function values at each of the n breakpoints (f). 
% The derivative is taken to be a piecewise linear function
% with values at the breakpoints equal to the 3-point finite
% difference approximations at the breakpoints. 
% These values are produced by D{end}*f.
%
% See also: LINBASE, LINDEF, FUNBASE, FUNDOP

% Copyright (c) 1997-2000, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function [D,n,a,b,parms]=lindop(breaks,evennum,order)

newbreaks=breaks;
n=length(breaks);
D=cell(abs(order),1);
for i=1:order
  d=1./diff(newbreaks);
  d=sparse([1:n-1 1:n-1],[1:n-1 2:n],[-d d],n-1,n);
  if i>1
    D{i}=d*D{i-1};
  else
    D{1}=d;
  end
  newbreaks=(newbreaks(1:end-1)+newbreaks(2:end))/2;
  n=n-1;
end
for i=-1:-1:order
  newbreaks=[[3 -1]*newbreaks(1:2);
    (newbreaks(1:end-1)+newbreaks(2:end));
    [-1 3]*newbreaks(end-1:end)]/2;
  d=diff(newbreaks)';
  n=n+1;
  d=tril(d(ones(n,1),:),-1);
  if i<-1
    D{-i}=d*D{-i-1};
  else
    D{1}=d;
  end
  % adjustment to make value at original left endpoint equal 0
  if evennum>0
    temp=linbase(newbreaks,length(newbreaks),breaks(1),0)*D{-i};
  else
    temp=linbase(newbreaks,0,breaks(1),0)*D{-i};
  end
  D{-i}=D{-i}-temp(ones(length(newbreaks),1),:);
end

if nargout>1
  n=length(newbreaks);
  a=newbreaks(1);
  b=newbreaks(end);
  if evennum>0
    parms={newbreaks,n};
  else
    parms={newbreaks,0};
  end
end