%% LCPLEMKE
%
%  Uses Lemke's algorithm to solve linear complementarity problem
%     z = M*x + q
%     a <= x <= b
%     x_i > a_i => z_i => 0
%     x_i < b_i => z_i =< 0
%
%  Usage
%    [x,z] = lcplemke(M,q,a,b,x)
%  Input
%    M       : n.n matrix
%    q       : n.1 vector
%    a       : n.1 vector, left bound on x
%    b       : n.1 vector, right bound on x
%    x       : n.1 vector, initial guess for solution
%  Output
%    x       : solution to lcp
%    z       : function value at x
%  Options
%    Options may be set with OPTSET
%    tol     : convergence tolerance
%    maxit   : maximum number of iterations

%  Copyright(c) 1997-2010
%    Mario J. Miranda - miranda.4@osu.edu
%    Paul L. Fackler  - paul_fackler@ncsu.edu

function [x,z] = lcplemke(M,q,a,b,x)

if nargin < 5, x=q; end

n    = length(x);
Minv = inv(M);
MM   = [-Minv Minv ; Minv -Minv];
qq   = [-Minv*q-a ; Minv*q+b];
z0   = [ max(-M*x-q,0) ; max(M*x+q,0) ];
z    = lemke(MM,qq,z0);
z    = z(n+1:2*n)-z(1:n);
x    = Minv*(z-q);


%% LEMKE 
%
%  Uses Lemke's algorithm to solve linear complementarity problem in
%  standard form
%     Mz+q >=0, z>=0, z'(Mz+q)=0.
%
%  Usage
%    [z,err] = lemke(M,q,z0)
%  Input
%    M    : n.n matrix
%    z    : n.1 vector
%    z0   : n.1 a starting basis (default: 0)
%  Note
%    The starting basis z0 can be either an initial guess of the solution 
%    or a vector of zeros and ones with ones representing those z(i) 
%    thought to be non-zero in the solution.  
%  Output
%    z    : n.1 solution
%    err  : 0-Solution found
%           1-Maximum iterations exceeded
%           2-Unbounded ray termination
%           3-Initial basis cannot be computed - try new value of z0
% Noete
%   Uses a modified Lemke's algorithm (complementary pivoting) modified to 
%   allow a user defined initial basis. 

function [z,err] = lemke(M,q,z0)

n = length(q);

if nargin<3, z0 = []; end

if size(M)~=[n n]
  error('Matrices are not compatible');
end

zer_tol = 1e-5;
piv_tol = 1e-8;
maxiter = min(1000,25*n);
err = 0;

% Trivial solution exists
if all(q>=0)
  z = zeros(n,1); return;
end

% Initializations 
z = zeros(2*n,1);
iter = 0;
x = q;

t = 2*n+1;                            % artificial variable

% Determine initial basis
if isempty(z0)
  bas = [];
  nonbas = (1:n)';
else
  bas = find(z0>0);
  nonbas = find(z0<=0);
end

B = spalloc(n,n,nnz(M));              % allocate memory for B
B = -speye(n);

% Determine initial values
if ~isempty(bas)
  B = [M(:,bas) B(:,nonbas)];
  if condest(B)>1e16
    z = []; 
    err = 3; 
    return
  end
  x = -(B\q);
end

% Check if initial basis provides solution
if all(x>=0)
  z(bas) = x(1:length(bas));
  z = z(1:n);
  return
end

% Determine initial leaving variable
[tval,lvindex] = max(-x);
bas = [bas;(n+nonbas)];
leaving = bas(lvindex);

bas(lvindex) = t;                     % pivot in the artificial variable

U = x<0;
Be = -(B*U);
x = x+tval*U;
x(lvindex) = tval;
B(:,lvindex) = Be;

% Main iterations begin here
for iter=1:maxiter
  % Check if done; if not, get new entering variable
  if (leaving == t)
    break
  elseif (leaving <= n)
    entering = n+leaving;
    Be = zeros(n,1); 
    Be(leaving) = -1;
  else
    entering = leaving - n;
    Be = M(:,entering);
  end
  d = B\Be;
  
  % Find new leaving variable
  j = find(d>piv_tol);                % indices of d>0
  if isempty(j)                       % no new pivots - ray termination
    err = 2;
    break
  end
  theta = min((x(j)+zer_tol)./d(j));  % minimal ratios, d>0
  j = j((x(j)./d(j))<=theta);         % indices of minimal ratios, d>0
  lvindex = find(bas(j)==t);          % check if artificial among these
  if ~isempty(lvindex)                % always use artifical if possible
    lvindex = j(lvindex);
  else                                % otherwise pick among set of max d
    theta = max(d(j));
    lvindex = find(d(j)==theta);
    lvindex = ceil(length(lvindex)*rand(1,1));  % if multiple choose randomly
    lvindex = j(lvindex);
  end
  leaving = bas(lvindex);
  
  % Perform pivot
  ratio = x(lvindex)./d(lvindex);
  x = x - ratio*d;
  x(lvindex) = ratio;
  B(:,lvindex) = Be;
  bas(lvindex) = entering;
end                                   % end of iterations
if iter>=maxiter && leaving~=t
  err = 1;
end

z(bas) = x; z = z(1:n);

% Display warning messages if no error code is returned
if nargout<2 && err~=0
  s = 'Warning: solution not found - ';
  if err==1
    disp([s 'Iterations exceeded limit']);
  elseif err==2
    disp([s 'Unbounded ray']);
  elseif err==3
    disp([s 'Initial basis infeasible']);
  end
end