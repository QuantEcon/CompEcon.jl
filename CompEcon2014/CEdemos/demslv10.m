function demslv10

%% DEMSLV10 Illustrates linear complementarity problem

% Preliminary tasks
demosetup(mfilename)

figure
basicsubplot(1,'f(a)>f(b)>0', 1.5, 0.5)
bullet(1.0, 0.5,18,'r')
basicsubplot(2,'f(a)>0>f(b)', 0.5,-0.5)
bullet(0.5, 0.0,18,'r')
basicsubplot(3,'0>f(a)>f(b)',-0.5,-1.5)
bullet(0.0,-0.5,18,'r')

figure
basicsubplot(1,'f(a)<f(b)<0',-1.5,-0.5)
bullet(0.0,-1.5,18,'r')
basicsubplot(2,'f(a)<0<f(b)',-0.5, 0.5)
bullet(0.0,-0.5,18,'r')
bullet(0.5, 0.0,18,'r')
bullet(1.0, 0.5,18,'r')
basicsubplot(3,'0<f(a)<f(b)', 0.5, 1.5)
bullet(1.0, 1.5,18,'r')

% Save Plots as EPS Files
printfigures(mfilename,2,1,0)


function basicsubplot(i,tit,yleft,yright)
fs = 14;  % FontSize
lw = 2;   % LineWidth
subplot(1,3,i)
hold on
title(tit,'Fontsize',fs)
plot([0;1],[yleft,yright],'LineWidth',5);
plot([0;1],[ 0;0],'k','LineWidth',lw);
plot([0;0],[-2;2],'--k','LineWidth',lw);
plot([1;1],[-2;2],'--k','LineWidth',lw);
ftext(-0.05,0,'0','right','middle',fs)
ftext(0,-2.8,'a','center','bottom',fs)
ftext(1,-2.8,'b','center','bottom',fs)
axis ([0 1 -2 2],'square','off')