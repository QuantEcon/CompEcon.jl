function demslv11

%% DEMSLV11 Minmax and semi-smooth rootfinding reformulations of
% complementarity problems.

% Preliminary tasks
demosetup(mfilename)
set(0,'DefaultLineLineWidth',2)

xmin = -0.3;
xmax =  0.3;
x = nodeunif(500,xmin,xmax);
z = zeros(size(x));

y = -4*x-2*tanh(x);
a = -0.5;
b =  0.5;
ya = a-x;
yb = b-x;

figure
for i=1:2
  subplot(1,2,i)
  hold on
  set(gca,'Xtick',[])
  set(gca,'Ytick',[])
  if i==1
    title('Minimax Reformulation')
    yr = minmax(x,a,b,y);
  else
    title('Semismooth Reformulation')
    yr = ssmooth(x,a,b,y);
  end
  plot(x,yr,'r','LineWidth',5);
  plot(x,y);
  plot(x,z,':');
  plot(x,ya,'k--');
  plot(x,yb,'k--');
  ftext(xmin,-0.3,'a-x' ,'left','top')
  ftext(xmax, 0.3,'b-x' ,'right','bottom')
  ftext(-0.15, 1.3,'f(x)','center','middle')
  axis ([xmin xmax -1.5 1.5],'square')
  box on
end

% Save Plots as EPS Files
printfigures(mfilename,1,1,0)