function demmath01

%% DEMMATH01 Taylor Approximations

% Preliminary tasks
demosetup(mfilename)

%% Univariate Taylor Approximation

x = nodeunif(100,-1,1);
y = (x+1).*exp(2*x);
y1 = 1+3*x;
y2 = 1+3*x+8*x.^2;

figure
hold on
plot(x,y ,'k','LineWidth',3)
plot(x,y1,'b','LineWidth',3)
plot(x,y2,'r','LineWidth',3)
legend('Function','1st Order Approximation','2nd Order Approximation','Location','nw')
legend 'boxoff'
set(gca,'XTick',[min(x) 0 max(x)])



%% Bivariate Taylor Approximation

nplot = [101 101];
a = [ 0 -1];
b = [ 2  1];
[x,xcoord] = nodeunif(nplot,a,b);

y  = exp(x(:,2)).*x(:,1).^2;
y1 = 2*x(:,1) - x(:,2) - 1;
y2 = x(:,1).^2 - 2*x(:,1).*x(:,2) + 0.5*x(:,2).^2 + x(:,2);
y  = reshape(y ,nplot(1),nplot(2));
y1 = reshape(y1,nplot(1),nplot(2));
y2 = reshape(y2,nplot(1),nplot(2));

figure
surf(xcoord{1},xcoord{2},y,'FaceColor','interp','EdgeColor','interp');
xlabel('x_1'); ylabel('x_2'); zlabel('y');
title('Original Function')

figure
surf(xcoord{1},xcoord{2},y1,'FaceColor','interp','EdgeColor','interp');
xlabel('x_1'); ylabel('x_2'); zlabel('y');
title('First-Order Approximation')

figure
surf(xcoord{1},xcoord{2},y2,'FaceColor','interp','EdgeColor','interp');
xlabel('x_1'); ylabel('x_2'); zlabel('y');
title('Second-Order Approximation')


%% Save Plots as EPS Files
printfigures(mfilename,4)