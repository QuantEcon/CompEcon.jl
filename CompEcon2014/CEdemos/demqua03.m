function demqua03

%% DEMQUA03  Area under 1-D and 2-D curves, various methods
%
% Uni- and bi-vaiariate integration using Newton-Cotes, Gaussian, Monte 
% Carlo, and quasi-Monte Carlo quadrature methods.

% Preliminary tasks
demosetup(mfilename)

figure
a = -1;
b =  1;
n = 301;
x = nodeunif(n,a,b);
y1 = exp(-x);
y2 = sqrt(abs(x));

subplot(1,2,1)
hold on
plot(x,y1)
set(gca,'xtick',[-1 0 1])
set(gca,'ytick',[0])
axis square
title('$e^{-x}$','interpreter','latex')

subplot(1,2,2)
hold on
plot(x,y2)
set(gca,'xtick',[-1 0 1])
set(gca,'ytick',[0])
axis square
title('$\sqrt{|x|}$','interpreter','latex')

% Save Plots as EPS Files
printfigures(mfilename,1)

n = [5;11;21;31];
N = length(n);
int = zeros(N,3);

% Define bounds
a = -1;
b =  1;

f = @(x) exp(-x);
for i=1:N
  [x1,w1] = qnwtrap(n(i),a,b);
  [x2,w2] = qnwsimp(n(i),a,b);
  [x3,w3] = qnwlege(n(i),a,b);
  int(i,1) = w1'*f(x1);
  int(i,2) = w2'*f(x2);
  int(i,3) = w3'*f(x3);
end
trueint = exp(1)-1/exp(1);
int = log10(abs(int/trueint-1));

% Print out 1-d integrations
fprintf('\n Integrate exp(-x) on [%1i,%1i]\n',a,b)
fprintf('         n    Trapezoid  Simpson    Legendre\n')
for i=1:length(n);
  fprintf('%10i %10.1f %10.1f %10.1f\n',[n(i) int(i,:)]);
end

f = @(x) abs(x).^0.5;
for i=1:N
  [x1,w1] = qnwtrap(n(i),a,b);
  [x2,w2] = qnwsimp(n(i),a,b);
  [x3,w3] = qnwlege(n(i),a,b);
  int(i,1) = w1'*f(x1);
  int(i,2) = w2'*f(x2);
  int(i,3) = w3'*f(x3);
end
trueint = 4/3;
int = log10(abs(int/trueint-1));

fprintf('\n Integrate abs(x).^0.5 on [%1i,%1i]\n',a,b)
fprintf('         n    Trapezoid  Simpson    Legendre\n')
for i=1:length(n);
  fprintf('%10i %10.1f %10.1f %10.1f\n',[n(i) int(i,:)]);
end