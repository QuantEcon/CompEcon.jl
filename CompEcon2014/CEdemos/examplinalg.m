%% Script file for linear algebra examples in lecture notes.


%% Dot product and angle

x = [0;1];
y = [1;1];
dotproduct = dot(x,y)
angle = acosd(dot(x,y)/sqrt(dot(x,x)*dot(y,y)))


%% Submatrices

A = [4 4;2 6;2 3];
A(3,:)
A(:,2)
A(3,1)


%% Matrix addition and subtraction

A = [4 4;2 6;2 3];
B = [1 6;5 2;3 4];
A+B
A-B


%% Matrix mutiplication
  
A = [1 5; 3 4; 6 3];
B = [3 1; 2 4];
A*B
% B*A


%% Matrix transposition

A = [2 2; 3 2; 1 0];
B = [1 2; 1 2; 3 1];
A''
(A+B)'
A'+B'
(A*B')'
B*A'


%% Nonstandard matrix operations

A = [1 5 4; 3 6 4];
B = [-1 1 2; 3 2 1];
A.*B
A./B
A.^B


%% Backslash operator

A = [1 1; 2 0];
b = [1; 2];
x = A\b
A*x


%% Matrix inverse 1

A = [1 3; 1 2];
inv(A)
inv(A)*A
A*inv(A)


%% Matrix inverse 2

A = [2 1; 1 0];
B = [1 1; 3 1];
inv(A*B)-inv(A)*inv(B)
inv(A*B)-inv(B)*inv(A)
inv(A')-inv(A)'


%% Determinants

A = [2 2 1; 3 1 2; 1 1 0];
B = [1 2 1; 2 3 1; 1 1 1];
det(A)
det(B)
det(A')
det(A*B)
det(3*A)
det(inv(A))


%% Eigenvalues and Eigenvectors 1

A = [4 5; -2 -3];
[X,D] = eig(A)
A = rand(5,5)
eig(A)
prod(eig(A))
det(A)
B = A+A'
eig(B)
prod(eig(B))
det(B)


%% Eigenvalues and Eigenvectors 2

A = [2 2 1; 3 1 2; 1 1 0];
A^3
[X,D] = eig(A);
X*(D^3)*inv(X)


%% Positive Definiteness and Cholesky Square Roots

A = [2 0 -1;0 3 4;-1 4 6];
eig(A)
U = chol(A)
U'*U