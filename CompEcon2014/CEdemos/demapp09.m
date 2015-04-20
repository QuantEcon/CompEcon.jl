function demapp09

%% DEMAPP09 Linear spline approximation

% Preliminary tasks
demosetup(mfilename)
fs = 12;

f = @(x) 50-cos(x.^2/8).*(x-pi+.5).^2;

xmin = 0;
xmax = 1.5*pi;
off = 0.05;
n = 401;
x = nodeunif(n,xmin,xmax);
y = f(x);
ymin = min(y);
ymax = max(y);
ywid = ymax-ymin;
ymin = ymin-0.5*ywid;
ymax = ymax+0.1*ywid;

figure
hold on
nnode = 4;
xnode = nodeunif(nnode,xmin,xmax);
ynode = f(xnode);
z = zeros(n,1);
for i=2:nnode
    j = find(x>=xnode(i-1)&x<=xnode(i));
    z(j) = ynode(i-1)+(x(j)-xnode(i-1))*(ynode(i)-ynode(i-1))/(xnode(i)-xnode(i-1));
end
plot(x,y)
plot(x,z,'r','linewidth',2)
xlim([xmin-off xmax+off])
ylim([ymin ymax])
for i=1:nnode
   plot([xnode(i) xnode(i)],[ymin ynode(i)],'k--','linewidth',2)
   ftext(xnode(i),ynode(i),'\bullet','center','middle',24)
end
set(gca,'xtick',[],'ytick',[])
ftext(xnode(1),ymin,'x_0=a','center','top',fs);
ftext(xnode(2),ymin,'x_1  ','center','top',fs);
ftext(xnode(3),ymin,'x_2  ','center','top',fs);
ftext(xnode(4),ymin,'x_3=b','center','top',fs);

figure
hold on
nnode = 7;
xnode = nodeunif(nnode,xmin,xmax);
ynode = f(xnode);
z = zeros(n,1);
for i=2:nnode
    j = find(x>=xnode(i-1)&x<=xnode(i));
    z(j) = ynode(i-1)+(x(j)-xnode(i-1))*(ynode(i)-ynode(i-1))/(xnode(i)-xnode(i-1));
end
plot(x,y)
plot(x,z,'r','linewidth',2)
xlim([xmin-off xmax+off])
ylim([ymin ymax])
for i=1:nnode
   plot([xnode(i) xnode(i)],[ymin ynode(i)],'k--','linewidth',2)
   ftext(xnode(i),ynode(i),'\bullet','center','middle',24)
end
set(gca,'xtick',[],'ytick',[])
ftext(xnode(1),ymin,'x_0=a','center','top',fs);
ftext(xnode(2),ymin,'x_1  ','center','top',fs);
ftext(xnode(3),ymin,'x_2  ','center','top',fs);
ftext(xnode(4),ymin,'x_3  ','center','top',fs);
ftext(xnode(5),ymin,'x_4  ','center','top',fs);
ftext(xnode(6),ymin,'x_5  ','center','top',fs);
ftext(xnode(7),ymin,'x_6=b','center','top',fs);

% Save Plots as EPS Files
printfigures(mfilename,2)