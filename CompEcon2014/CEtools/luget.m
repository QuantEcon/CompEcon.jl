% luget LU factorization for use with lusol
% USAGE
%   ALU=luget(A,alg);
% INPUTS
%   A   : nxn matrix
%   alg : optional flag for sparse A; alg controls algorithm options
%           2 : only L and U are computed (L is not triangular)
%           3 : L, U and a row permutation vector are computed
%           4 : columns of A are permuted before factorization
%           5 : A scaled and column permuted before factorization [default]
% OUTPUT
%   ALU : cell array containing LU factorization
%
% ALU{1}=L, ALU{2}=U, ALU{3}=P or =P/R, ALU{4}=Q
% For a description of these matrices see lu
%
% To solve Ax=b use
%   ALU=luget(A); x=lusol(ALU,b);
% Note that this is basically equivalent to x=A\b except that the latter
% refines x using iterative refinement. To accomplish this one could use
%   ALU=luget(A); x=lusol(ALU,b); x=x+lusol(ALU,b-A*x);  
%
% The algorithms used by lu (and thus by luget) can be controled using spparms
%
% See: lusol, lu

function ALU=luget(A,alg)

n=size(A);
if ndims(n)>2 || n(1)~=n(2)
  error('A must be a square matrix');
end

if issparse(A)
  if nargin<2 || isempty(alg), alg=5; end
  switch alg
    case 5
      ALU=cell(1,5);
      [ALU{:}]=lu(A);
      % combine row permutation and row scaling matrices
      ALU{3}=ALU{3}/ALU{5}; ALU(5)=[]; 
    case 4
      ALU=cell(1,4);
      [ALU{:}]=lu(A);
    case 3
      ALU=cell(1,3);
      [ALU{:}]=lu(A);
    case 2
      ALU=cell(1,2);
      [ALU{:}]=lu(A);
    otherwise
      error('Incorrect specification of alg')
  end
else
  ALU=cell(1,3);
  [ALU{:}]=lu(A);
end