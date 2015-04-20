function demopt01

%% DEMOPT01 Maximization via golden search

% Preliminary tasks
demosetup(mfilename)

x = nodeunif(500,0,3);
f = @(x) x.*cos(x.^2);

x1 = golden(f,0,1);
x2 = golden(f,2,3);

figure
hold on
plot(x,f(x),'k')
bullet(x1,f(x1),25,'r')
bullet(x2,f(x2),25,'g')
title('Maximization of x cos(x^2) via golden search')

% Save Plots as EPS Files
printfigures(mfilename,1)