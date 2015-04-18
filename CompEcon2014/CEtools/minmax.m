% MINMAX  Minimax transformation for solving NCP as rootfinding problem
% USAGE
%   [fhatval,fhatjac] = minmax(x,a,b,fval,fjac)
% INPUTS
%   x       : n by 1 vector, evaluation point
%   a       : n by 1 vector, left bound on x
%   b       : n by 1 vector, right bound on x
%   fval    : function value at x
%   fjac    : function Jacobian at x
% OUTPUTS
%   fhatval : transformed function value at x
%   fhatjac : transformed function Jacobian at x
%
% See also: smooth

% Copyright (c) 1997-2000, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function [fhatval,fhatjac] = minmax(x,a,b,fval,fjac)

if length(a)==1,a=a*ones(size(x)); end
if length(b)==1,b=b*ones(size(x)); end

da = a-x;
db = b-x;
fhatval = min(max(fval,da),db); 

if nargout==2
   fhatjac = -eye(length(x));
   i = find(fval>da & fval<db);
   fhatjac(i,:) = fjac(i,:);
end
