% FUNFITF  Computes interpolation coefficients for D-dim function.
% Approximates a function y=f(x) using a specified basis type.
% USAGE
%   c=FUNFITF(basis,f,varargin);
% INPUTS
%   basis    : a structure defining a family of functions (use FUNDEF to create this)
%   f         : a function with a single input argument, x and single output argument y.
%   [varargin]: additional parameters passed to f
%  OUTPUTS
%    c   : a coefficient vector or matrix
%    B   : the basis structure used to evaluate the function
%
% The approximating function can then be evaluated using FUNEVAL(G,X)
%
% f may be either a an 'inline' variable or the name of a function file.
%   In either case it should be accept a kxd matrix and
%      kxp matrix of function values at the x(i,:), i=1..k.
%   Thus f is a mapping from a d-dimensional to a p-dimensional space.
%
% Example:
%   f=inline(exp(x(:,1).*x(:,2)),'x');
%   basis=fundef({'cheb',10,0,1},{'cheb',6,-1,1});
%   c=funfitf(basis,f);
% Alternatively one can define a function file F.m:
%   function y=F(x)
%   y=exp(x(:,1).*x(:,2));
% and call FUNFITF using
%   c=funfitf(basis,'F');
% (note that F is enclosed in quotes).
% Either method yields a coefficient structure that can be
%    evaluated using FUNEVAL.
%
% To fit a function to data use FUNFITXY.
%
% USES: funnode, funbasex, ckronxi
%
% See also: FUNFITXY, FUNDEF, FUNEVAL, FUNNODE, FUNBASE.

% Copyright (c) 1997-2000, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function c=funfitf(basis,f,varargin)

if nargin<2 error('At least two parameters must be specified'); end
if nargin<3 varargin={}; end

x = funnode(basis);
y = feval(f,x,varargin{:});
c = funfitxy(basis,x,y);