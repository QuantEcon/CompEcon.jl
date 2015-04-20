function demqua08

%% DEMQUA08 Illustrates integration using Simpson's rule

% Preliminary tasks
demosetup(mfilename)


%% Function Definition
c = [0.50; -1.00; 2.00];
f = @(x) [x x.^2 x.^3]*c;


%% Basic Figure Setup
xmin = -1;
xmax =  1;
xwid = xmax-xmin;
n = 401;
x = nodeunif(n,xmin,xmax);
y = f(x);
ymin = min(y);
ymax = max(y);
ywid = ymax-ymin;
ymin = ymin-0.2*ywid;
ymax = ymax+0.1*ywid;


%% Two Intervals

figure
hold on
nnode = 3;
xnode = nodeunif(nnode,xmin,xmax);
c = [xnode.^0 xnode.^1 xnode.^2]\f(xnode);
g = @(x) [x.^0 x.^1 x.^2]*c;
y = g(x);
for i=1:n
   plot([x(i) x(i)],[ymin+0.02 y(i)],'y-','linewidth',2)
end
plot(x,f(x))
plot(x,y,'r--','linewidth',2)
xlim([xmin-0.1*xwid xmax+0.1*xwid])
ylim([ymin ymax])
for i=1:nnode
   plot([xnode(i) xnode(i)],[ymin f(xnode(i))],':','linewidth',1)
   ftext(xnode(i),f(xnode(i)),'\bullet','center','middle');
end
set(gca,'xtick',[])
set(gca,'ytick',[])
ftext(xnode(1),ymin-0.1,'x_0=a','center','top');
ftext(xnode(2),ymin-0.1,'x_1  ','center','top');
ftext(xnode(3),ymin-0.1,'x_2=b','center','top');


%% Four Intervals
figure
hold on

m = (n+1)/2;
x1 = x(1:m);
x2 = x(m:n);
nnode = 5;
mnode = (nnode+1)/2;
xnode = nodeunif(nnode,xmin,xmax);
xnode1 = xnode(1    :mnode);
xnode2 = xnode(mnode:nnode);
c1 = [xnode1.^0 xnode1.^1 xnode1.^2]\f(xnode1);
c2 = [xnode2.^0 xnode2.^1 xnode2.^2]\f(xnode2);
g1 = @(x) [x.^0 x.^1 x.^2]*c1;
g2 = @(x) [x.^0 x.^1 x.^2]*c2;
y = [g1(x1);g2(x2)]; y(m) = [];
for i=1:n
   plot([x(i) x(i)],[ymin+0.02 y(i)],'y-','linewidth',2)
end
plot(x,f(x))
plot(x,y,'r--','linewidth',2)
xlim([xmin-0.1*xwid xmax+0.1*xwid])
ylim([ymin ymax])
for i=1:nnode
   plot([xnode(i) xnode(i)],[ymin f(xnode(i))],':','linewidth',1)
   ftext(xnode(i),f(xnode(i)),'\bullet','center','middle');
end
set(gca,'xtick',[])
set(gca,'ytick',[])
ftext(xnode(1),ymin-0.1,'x_0=a','center','top');
ftext(xnode(2),ymin-0.1,'x_1  ','center','top');
ftext(xnode(3),ymin-0.1,'x_2  ','center','top');
ftext(xnode(4),ymin-0.1,'x_3  ','center','top');
ftext(xnode(5),ymin-0.1,'x_4=b','center','top');

%% Save Plots as EPS Files
printfigures(mfilename,2)