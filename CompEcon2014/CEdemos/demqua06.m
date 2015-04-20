function demqua06

%% DEMQUA06 % Area under a curve

% Preliminary tasks
demosetup(mfilename)


figure
hold on
xmin = 0;
xmax = 1;
a = 0.25;
b = 0.75;
n = 401;
ymin = 25;
ymax = 65;
x = nodeunif(n,a,b);
y = f(x);
for i=1:n
   plot([x(i) x(i)],[ymin+0.015 y(i)],'y-','linewidth',2)
end
plot([a a],[ymin f(a)],'k--','linewidth',2)
plot([b b],[ymin f(b)],'k--','linewidth',2)
x = nodeunif(n,xmin,xmax);
y = f(x);
plot(x,y)
ylim([ymin ymax])
set(gca,'xtick',[])
set(gca,'ytick',[])
ftext(xmin,ymin,'0','right','top');
ftext(a,ymin,'a','center','top');
ftext(b,ymin,'b','center','top');
ftext(xmax,f(xmax)-0.2,'f','center','top');
ftext((xmin+xmax)/2,(ymin+ymax)/2.5,'\int_a^b f(x)dx','center','middle');

% Save Plots as EPS Files
printfigures(mfilename,1)


function y=f(x)
y = 50-cos(pi*x).*(2*pi*x-pi+.5).^2;