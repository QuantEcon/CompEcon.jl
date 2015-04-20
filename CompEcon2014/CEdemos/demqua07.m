function demqua07

%% DEMQUA07 Illustrates integration using Trapezoidal rule

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
ynode = f(xnode);
z = zeros(n,1);
for i=2:nnode
    j = find(x>=xnode(i-1)&x<=xnode(i));
    z(j) = ynode(i-1)+(x(j)-xnode(i-1))*(ynode(i)-ynode(i-1))/(xnode(i)-xnode(i-1));
end
for i=1:n
   plot([x(i) x(i)],[ymin+0.02 z(i)],'y-','linewidth',2)
end
plot(x,y)
plot(x,z,'r--','linewidth',2)
xlim([xmin-0.1*xwid xmax+0.1*xwid])
ylim([ymin ymax])
for i=1:nnode
   plot([xnode(i) xnode(i)],[ymin ynode(i)],':','linewidth',1)
   ftext(xnode(i),ynode(i),'\bullet','center','middle');
end
set(gca,'xtick',[])
set(gca,'ytick',[])
ftext(xnode(1),ymin-0.1,'x_0=a','center','top');
ftext(xnode(2),ymin-0.1,'x_1  ','center','top');
ftext(xnode(3),ymin-0.1,'x_2=b','center','top');


%% Four Intervals
figure
hold on
nnode = 5;
xnode = nodeunif(nnode,xmin,xmax);
ynode = f(xnode);
z = zeros(n,1);
for i=2:nnode
    j = find(x>=xnode(i-1)&x<=xnode(i));
    z(j) = ynode(i-1)+(x(j)-xnode(i-1))*(ynode(i)-ynode(i-1))/(xnode(i)-xnode(i-1));
end
for i=1:n
   plot([x(i) x(i)],[ymin+0.02 z(i)],'y-','linewidth',2)
end
plot(x,y)
plot(x,z,'r--','linewidth',2)
xlim([xmin-0.1*xwid xmax+0.1*xwid])
ylim([ymin ymax])
for i=1:nnode
   plot([xnode(i) xnode(i)],[ymin ynode(i)],':','linewidth',1)
   ftext(xnode(i),ynode(i),'\bullet','center','middle');
end
set(gca,'xtick',[])
set(gca,'ytick',[])
ftext(xnode(1),ymin-0.1,'x_0=a','center','top');
ftext(xnode(2),ymin-0.1,'x_1  ','center','top');
ftext(xnode(3),ymin-0.1,'x_2  ','center','top');
ftext(xnode(4),ymin-0.1,'x_3  ','center','top');
ftext(xnode(5),ymin-0.1,'x_4=b','center','top');


%% Eight Intervals
figure
hold on
nnode = 9;
xnode = nodeunif(nnode,xmin,xmax);
ynode = f(xnode);
z = zeros(n,1);
for i=2:nnode
    j = find(x>=xnode(i-1)&x<=xnode(i));
    z(j) = ynode(i-1)+(x(j)-xnode(i-1))*(ynode(i)-ynode(i-1))/(xnode(i)-xnode(i-1));
end
for i=1:n
   plot([x(i) x(i)],[ymin+0.02 z(i)],'y-','linewidth',2)
end
plot(x,y)
plot(x,z,'r--','linewidth',2)
xlim([xmin-0.1*xwid xmax+0.1*xwid])
ylim([ymin ymax])
for i=1:nnode
   plot([xnode(i) xnode(i)],[ymin ynode(i)],':','linewidth',1)
   ftext(xnode(i),ynode(i),'\bullet','center','middle');
end
set(gca,'xtick',[])
set(gca,'ytick',[])
ftext(xnode(1),ymin-0.1,'x_0=a','center','top');
ftext(xnode(2),ymin-0.1,'x_1  ','center','top');
ftext(xnode(3),ymin-0.1,'x_2  ','center','top');
ftext(xnode(4),ymin-0.1,'x_3  ','center','top');
ftext(xnode(5),ymin-0.1,'x_4  ','center','top');
ftext(xnode(6),ymin-0.1,'x_5  ','center','top');
ftext(xnode(7),ymin-0.1,'x_6  ','center','top');
ftext(xnode(8),ymin-0.1,'x_7  ','center','top');
ftext(xnode(9),ymin-0.1,'x_8=b','center','top');


%% Save Plots as EPS Files
printfigures(mfilename,3)