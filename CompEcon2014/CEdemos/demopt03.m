function demopt03

%% DEMOPT03 Nelder-Mead simplex method
%
% Creates and plays a movie of Nelder-Meade simplex iterations when
% maximizing banana function f(x,y)=-100*(y-x*x)^2-(1-x)^2, starting at [0;1].

% Preliminary tasks
demosetup(mfilename)
warning off

%  The "banana" or Rosencrantz function
banana = @(x) -100*(x(2,:)-x(1,:).^2).^2-(1-x(1,:)).^2;

n = [20 20];
xmin = [-0.2 -0.2];
xmax = [ 1.2  1.2];
[x,xcoord] = nodeunif(n,xmin,xmax);

y = banana(x');
y = reshape(y,n)';
conts = -exp(0.25:0.5:20);

figure
hold on
contour(xcoord{1},xcoord{2},y,conts,'k:')
xlabel('x_1')
ylabel('x_2')
title('Nelder-Mead Maximizes the Banana Function')
set(gca,'ytick',[0 0.5 1])
axis([-.2 1.2 -.2 1.2])

optset('neldmead','maxit',1);
x = [1;0];
[xx,S] = neldmead(banana,x);
hp = patch(S(1,:),S(2,:),[0.5 0.5 0.5]);
for i=1:60
  xvec(:,i) = x;
  [x,S] = neldmead(banana,x,S);
  set(hp,'xdata',S(1,:)','ydata',S(2,:)');
  getframe;
end