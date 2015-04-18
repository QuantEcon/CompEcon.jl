% LUSOL Sparse LU solver for Ax=b
% USAGE
%   x=lusol(ALU,b);
% INPUTS
%   ALU      : cell array containing factorization of A
%                obtained from luget
%   b        : dense matrix (n x m)
% OUTPUT
%   x        : dense matrix (n x m)
%
% To solve Ax=b use
%   ALU=luget(A); x=lusol(ALU,b);
% Note that this is basically equivalent to x=A\b except that the latter
% refines x using iterative refinement. To accomplish this one could use
%   ALU=luget(A); x=lusol(ALU,b); x=x+lusol(ALU,b-A*x);  
%
% see luget for description of ALU
function x=lusol(ALU,b)
if ~isa(ALU,'cell')
  error('ALU must be a cell array: use luget to define ALU')
end
m=length(ALU);
switch m
      case 2
        x=ALU{2}\(ALU{1}\b);
      case 3
        x=ALU{2}\(ALU{1}\(ALU{3}*b));
      case 4
        x=ALU{4}*(ALU{2}\(ALU{1}\(ALU{3}*b)));
      otherwise
        error('LU information is not properly specified')
end