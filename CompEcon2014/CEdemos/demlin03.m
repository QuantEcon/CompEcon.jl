function demlin03

%% DEMLIN03 - Sparse linear equations


% Preliminary tasks
demosetup(mfilename)

AA = rand(1000,1000);
bb = rand(1000,1);
for i=1:1000
  for j=1:1000
    if abs(i-j)>1, AA(i,j) = 0; end
  end
end

n = [50:50:200 300:100:1000];
for k=1:length(n)
  A = AA(1:n(k),1:n(k));
  b = bb(1:n(k));
  tic
  for i=1:100
    x = A\b;
  end
  toc1 = toc;
  S = sparse(A);
  tic
  for i=1:100
    x = S\b;
  end
  toc2 = toc;
  ratio(k) = toc2/toc1;
end

% Plot effort ratio
figure
plot(n,ratio)
xlabel('n'); 
ylabel('Ratio');
title('Ratio of Sparse Solve Time to Full Solve Time')

% Save Plots as EPS Files
printfigures(mfilename,1)