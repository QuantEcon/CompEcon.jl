function demopt06

%% DEMOPT06 KKT conditions for constrained optimization problems

% Preliminary tasks
demosetup(mfilename)

x=nodeunif(100,-0.5,1.5);
a=0.1;
b=0.9;

figure

subplot(1,2,1)
f = @(x) 1.5-2*(x-.75).^2;
hold on;
plot(x,f(x));
plot([a;a],[-0.5;2],'k','LineWidth',2);
plot([b;b],[-0.5;2],'k','LineWidth',2);
set(gca,'XTick',[a b])
set(gca,'XTickLabel',{'a' 'b'})
set(gca,'YTick',[])
xstar = 0.75;
ystar = f(xstar);
bullet(xstar,ystar)
axis([a b 0.5 2.0],'square')
title('Internal Maximum')

subplot(1,2,2)
f = @(x) 1-(x+0.25).^2;
hold on;
plot(x,f(x));
plot([a;a],[-0.5;2],'k','LineWidth',2);
plot([b;b],[-0.5;2],'k','LineWidth',2);
set(gca,'XTick',[a b])
set(gca,'XTickLabel',{'a' 'b'})
set(gca,'YTick',[])
xstar = a;
ystar = f(xstar);
bullet(xstar,ystar)
axis([a b -0.5 1.0],'square')
title('Corner Maximum')

% Save Plots as EPS Files
printfigures(mfilename,1,1,0)