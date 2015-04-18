% FUNFITXY  Computes interpolation coefficients for d-dim function.
% USAGE
%   [c,B] = funfitxy(basis,x,y);
% or
%   [c,B] = funfitxy(basis,B,y);
% Approximates a function y = f(x) using a specified basis type.
% INPUTS
%   basis  : a structure defining a family of functions (use FUNDEF to create this)
%   x       : input values (kxd matrix or 1xd cell array of vectors)
%   B       : a basis structure - if funfitxy is called multiple time
%               the basis information can be stored and reused (see below)
%   y       : kxp matrix of output values
% OUTPUTS
%   c   : a coefficient matrix
%   B   : the basis structure or matrix used to evaluate the function
%
% The approximating function can then be evaluated using 
%   funeval(g,basis,X)
%
% Note: when evaluating on a multidimensional grid, pass x as a 
% 1 x d cell array of values. An efficient algorithm can be 
% used in this case.
%
% Example:
%   basis = fundefx({'cheb',10,0,1},{'cheb',6,-1,1});
%   g = funfitxy(basis,x,y);
% This produces a 2-d polynomial approximation to the data (x,y),
% where x must be a 2-column matrix or a 1x2 cell array.
% If the number of evaluation points is equal to prod(N) (60 in the
% example above) the approximation will interpolate those points.
% If there are more than this number of evaluation points, a least
% squares fit is provided.
%
% Reusing the basis information:
%   [c1,B] = funfitxy(basis,x,y1);
%   c2 = funfitxy(basis,B,y2);
% produces more efficiently the same result as
%   c1 = funfitxy(basis,x,y1);
%   c2 = funfitxy(basis,x,y2);
% Note: B.order(1,:) must equal 0.
%
% To fit a function using a MATLAB function (M-file or inline function) use FUNFITF.
%
% USES: FUNBASEX, CKRONXI
%
% See also: FUNFITF, FUNDEF, FUNEVAL, FUNNODE, FUNBASE.

% Copyright (c) 1997-2000, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function [c,B] = funfitxy(basis,x,y)

if nargin~=3, error('Three parameters must be specified'); end

m = size(y,1);
% There cannot be more basis functions than data points
if prod(basis.n)>m
    error('Insufficient number of data points')
end

if isstruct(x)                 % Use precomputed basis structure
    B = x;
    switch B.format
        case 'tensor'
            if any(B.order(1,:)~=0)
                error('invalid basis structure - first elements must have order 0')
            end
            B.vals = B.vals(1,:);
            c = ckronxi(B.vals,y,basis.d:-1:1);
        case 'direct'
            B = funbconv(B,zeros(1,size(B.vals,2)),'expanded');
            c = B.vals{1}\y;
        case 'expanded'
            c = B.vals{1}\y;
    end
elseif iscell(x)               % evaluate at grid points
    mm = 1;for i=1:size(x,2), mm=mm*size(x{i},1); end
    if mm~=m
        error('In FUNFITXY: x and y are incompatible')
    end
    B = funbasex(basis,x,0);
    c = ckronxi(B.vals,y,basis.d:-1:1);
else                           % evaluate at arbitrary points
    % x and y must have the same number of rows
    if size(x,1)~=m
        error('In FUNFITXY: x and y are incompatible')
    end
    B = funbasex(basis,x,0,'expanded');
    c = B.vals{1}\y;
end  