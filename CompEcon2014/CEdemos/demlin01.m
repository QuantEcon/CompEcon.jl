function demlin01

%% DEMLIN01 Solving linear equations by different methods

% Preliminary tasks
demosetup(mfilename)

% Print table header
fprintf('Hundreds of seconds required to solve n by n linear equation Ax=b\n')
fprintf('m times using A\\b and inv(A)*b, computing inverse only once.\n')
fprintf('    m       n        A\\b   inv(A)*b\n')

for m=[1,100]
for n=[50,500]
   A = rand(n,n);
   b = rand(n,1);
   tic
   for j=1:m
      x = A\b;
   end
   f1 = 100*toc;
   tic
   Ainv = inv(A);
   for j=1:m
      x = Ainv*b;
   end
   f2 = 100*toc;
   fprintf('  %3i     %3i  %9.2f  %9.2f \n',m,n,f1,f2)
end
end