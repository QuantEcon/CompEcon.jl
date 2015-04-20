function demqua01

%% DEMQUA01 Equidistributed sequences on unit square in R^2

% Preliminary tasks
demosetup(mfilename)

% Number of integration nodes
n = 2500;

% Neiderreiter Sequence
figure
x = qnwequi(n,[0 0],[1 1],'N');
plot(x(:,1),x(:,2),'.','MarkerSize',7);
title('Neiderreiter Sequence')
xlabel('x_1')
ylabel('x_2')
axis square
set(gca,'xtick',[0 1],'ytick',[0 1])

% Weyl Sequence
figure
x = qnwequi(n,[0 0],[1 1],'W');
plot(x(:,1),x(:,2),'.','MarkerSize',7);
title('Weyl Sequence')
xlabel('x_1')
ylabel('x_2')
axis square
set(gca,'xtick',[0 1],'ytick',[0 1])

% Pseudo-Random Sequence
figure
x = qnwequi(n,[0 0],[1 1],'R');
plot(x(:,1),x(:,2),'.','MarkerSize',7);
title('Pseudo-Random Sequence')
xlabel('x_1')
ylabel('x_2')
axis square
set(gca,'xtick',[0 1],'ytick',[0 1])

% Save Plots as EPS Files
printfigures(mfilename,3)