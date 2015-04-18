%% CHECKJAC 
%
%  Compares analytic and finite difference derivatives
%
%  Usage
%    [error,i,j] = checkjac(f,x,varargin)
%  Input
%    f        : name of function of form [fval,fjac]=f(x)
%               where fval and fjac are value and jacobian of f
%    x        : evaluation point
%    varargin : parameters passed to function f
%  Output
%    error    : maximum difference between analytic and finite difference Jacobians
%    i,j      : indices of maximum difference
% 
%  Uses: fdjac

%  Copyright(c) 1997-2010
%   Mario J. Miranda - miranda.4@osu.edu
%   Paul L. Fackler  - paul_fackler@ncsu.edu

function [error,i,j] = checkjac(f,x,varargin)

[fval,fjac] = feval(f,x,varargin{:});
fjacfindif  = fdjac(f,x,varargin{:});
[error,i]   = max(abs(fjac-fjacfindif));
[error,j]   = max(error);
i = i(j);
fprintf('Checking Jacobian\n')
fprintf('  Max Error  %6.2g\n',error)
fprintf('  Row        %7i\n',i)
fprintf('  Column     %7i\n',j)