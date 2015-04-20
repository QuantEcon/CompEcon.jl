function demslv06

%% DEMSLV06 Illustrates function iteration, Newton, and secant methods

% Preliminary tasks
demosetup(mfilename)


%% Function Iteration

xmin = 0.0;
xmax = 1.4;
xinit = 0.3;
xstar = 0.5*(1+sqrt(1.8));
xx = linspace(xmin,xmax);
yy = g(xx);

figure
hold on
plot(xx,yy)
plot(xx,xx,'k-','LineWidth',2)
axis([min(xx) max(xx) min(xx) max(xx)])
title('Function Iteration')

x(1) = xinit;
y(1) = g(x(1));
for i=2:20
  x(i) = y(i-1);
  y(i) = g(x(i));
end
ftext(x(1),0,'x_0','center','top')
ftext(x(2),0,'x_1','center','top')
ftext(x(3),0,'x_2','center','top')
ftext(xstar,0,'x*','center','top')
ftext(0,x(2),'x_1','right','middle')
ftext(0,x(3),'x_2','right','middle')
ftext(0,x(4),'x_3','right','middle')
ftext(0,xstar,'x*','right','middle')
for i=[1:3 19]
  plot([x(i) x(i)],[0 x(i+1)],'k--','LineWidth',2);
  plot([0 x(i)],[x(i+1) x(i+1)],'k--','LineWidth',2);
end
for i=2:20
  plot([x(i-1) x(i-1)],[x(i-1) x(i)],'r','LineWidth',3)
  plot([x(i-1) x(i)],[x(i) x(i)],'r','LineWidth',3)
end
for i=1:3
  bullet(x(i),x(i),18)
  bullet(x(i),x(i+1),18)
end
ftext(0.15,0.01,'45^o','center','bottom',12)
ftext(xmax,sqrt(xmax+.2),'g','center','top')
bullet(xstar,xstar,18)
set(gca,'XTick',[])
set(gca,'YTick',[])


%% Newton's Method

xmin = 1.0;
xmax = 2.5;
xinit = xmax;
xstar = 3^(1/5);
xx = linspace(xmin,xmax);
yy = f(xx);
zz = zeros(size(xx));

figure
hold on
plot(xx,yy)
plot(xx,zz,'k')
title('Newton''s Method')
axis off

x(1) = xinit;
y(1) = f(x(1));
for i=2:4
  [ylag dlag] = f(x(i-1));
  x(i) = x(i-1) - ylag/dlag;
  y(i) = f(x(i));
end
ftext(x(1),0,'x_0','center','top')
ftext(x(2),0,'x_1','center','top')
ftext(x(3),0,'x_2','center','top')
ftext(x(4),0,'x_3','center','top')
ftext(xstar,0,'x*','center','top')
for i=2:4
  plot([x(i-1) x(i)],[y(i-1) 0],'r','LineWidth',3)
end
for i=1:4
  plot([x(i) x(i)],[0 y(i)],'k--','LineWidth',2)
  bullet(x(i),y(i),18)
  bullet(x(i),0,18)
end
bullet(xstar,0,18)
ftext(xmin-.02,0,'0','right','middle')
ftext(xmax,f(xmax),'f','center','bottom')
set(gca,'XTick',[])
set(gca,'YTick',[])


%% Secant Method

figure
hold on
plot(xx,yy);
plot(xx,zz,'k');
title('Secant Method')
axis off

x(1) = xinit;
x(2) = xinit-0.25;
y(1) = f(x(1));
y(2) = f(x(2));
for i=3:4
  x(i) = x(i-1)-y(i-1)*(x(i-1)-x(i-2))/(y(i-1)-y(i-2));
  y(i) = f(x(i));
end
ftext(x(1),0,'x_0','center','top')
ftext(x(2),0,'x_1','center','top')
ftext(x(3),0,'x_2','center','top')
ftext(x(4),0,'x_3','center','top')
ftext(xstar,0,'x*','center','top')
for i=3:4
  plot([x(i) x(i-2)],[0 y(i-2)],'r','LineWidth',3);
end
for i=1:4
  plot([x(i) x(i)],[0 y(i)],'k--','LineWidth',2)
  bullet(x(i),y(i),18)
  bullet(x(i),0,18)
end
bullet(xstar,0,18)
ftext(xmin-.02,0,'0','right','middle')
ftext(xmax,f(xmax),'f','center','bottom')
set(gca,'XTick',[])
set(gca,'YTick',[])

% Save Plots as EPS Files
printfigures(mfilename,3)


%% Function
function [fval,fjac] = f(x)
fval = x.^5-3;
fjac = 5*x.^4;

%% Function
function gval = g(x)
gval = (x+0.2).^0.5';