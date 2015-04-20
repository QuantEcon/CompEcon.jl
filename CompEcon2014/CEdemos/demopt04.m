function demopt04

%% DEMOPT04 Maximization of banana function by various methods
% f(x,y)=-100*(y-x*x)^2-(1-x)^2, starting at [0;1].

% Preliminary tasks
demosetup(mfilename)
warning off

optset('qnewton','maxit',250);
optset('qnewton','ShowIters',0);

%  The "banana" or Rosencrantz function
banana  = @(x) -100*(x(2,:)-x(1,:).^2).^2-(1-x(1,:)).^2;

% for i=1:3
%     for j=1:4
%         optset('qnewton','SearchMeth',i);
%         optset('qnewton','StepMeth',j);
%         x0 = [1;0];
%         x = qnewton(banana,x0);
%         err = norm(x-[1;1]);
%         fprintf('%4i %4i  %12.8f\n',[i;j;err])
%     end
% end

optset('qnewton','ShowIters',1);
optset('neldmead','ShowIters',1);
optset('neldmead','maxit',250);

'Steepest Ascent Maximization'
optset('qnewton','SearchMeth',1);
x0 = [1;0];
x = qnewton(banana,x0)

'DFP Maximization'
optset('qnewton','SearchMeth',2);
x0 = [1;0];
x = qnewton(banana,x0)

'BFGS Maximization'
optset('qnewton','SearchMeth',3);
x0 = [1;0];
x = qnewton(banana,x0)

J = fdjac(banana,x)
H = fdhess(banana,x)
E = eig(H)

'Nelder-Mead Maximization'
x0 = [1;0];
x = neldmead(banana,x0)

optset('qnewton','maxit',1);
optset('neldmead','maxit',1);
optset('qnewton','ShowIters',0);
optset('neldmead','ShowIters',0);

n = [40 40];
xmin = [-0.7 -0.2];
xmax = [ 1.2  1.2];
[x,xcoord] = nodeunif(n,xmin,xmax);
y = banana(x');
y = reshape(y,n(1),n(2))';
conts = -exp(0.25:0.5:20);


%% Steepest Ascent Maximization

optset('qnewton','SearchMeth',1);
k = 250;
x = [1;0];
[xtmp,A] = qnewton(banana,x);
for i=1:k
  xx1(:,i) = x;
  [x,A] = qnewton(banana,x,A);
  if norm(x-[1;1])>sqrt(eps), iters1=i; end
end

figure
hold on
plot(xx1(1,:),xx1(2,:))
bullet(xx1(1,:),xx1(2,:),15,'r')
bullet(1,1,25,'k')
contour(xcoord{1},xcoord{2},y,conts,'k:')
xlabel('x_1','VerticalAlignment','cap')
ylabel('x_2','VerticalAlignment','bottom')
title('Steepest Ascent Maximization of Banana Function')
set(gca,'ytick',[0 0.5 1])
axis([-.7 1.2 -.2 1.2])


%% DFP Maximization

optset('qnewton','SearchMeth',2);
x = [1;0];
[xtmp,A] = qnewton(banana,x);
for i=1:k
  xx2(:,i) = x;
  [x,A] = qnewton(banana,x,A);
  if norm(x-[1;1])>sqrt(eps), iters2=i; end
end

figure
hold on
plot(xx2(1,:),xx2(2,:))
bullet(xx2(1,:),xx2(2,:),15,'r')
bullet(1,1,25,'k')
contour(xcoord{1},xcoord{2},y,conts,'k:')
xlabel('x_1')
ylabel('x_2')
title('DFP Quasi-Newton Maximization of Banana Function')
set(gca,'ytick',[0 0.5 1])
axis([-.7 1.2 -.2 1.2])


%% BFGS Maximization

optset('qnewton','SearchMeth',3);
x = [1;0];
[xtmp,A] = qnewton(banana,x);
for i=1:k
  xx3(:,i) = x;
  [x,A] = qnewton(banana,x,A);
  if norm(x-[1;1])>sqrt(eps), iters3=i; end
end

figure
hold on
plot(xx3(1,:),xx3(2,:))
bullet(xx3(1,:),xx3(2,:),15,'r')
bullet(1,1,25,'k')
contour(xcoord{1},xcoord{2},y,conts,'k:')
xlabel('x_1')
ylabel('x_2')
title('BFGS Quasi-Newton Maximization of Banana Function')
set(gca,'ytick',[0 0.5 1])
axis([-.7 1.2 -.2 1.2])


%% Nelder-Mead Maximization

x = [1;0];
[xx,S] = neldmead(banana,x);
for i=1:60
  xx4(:,i) = x;
  [x,S] = neldmead(banana,x,S);
  if norm(x-[1;1])>sqrt(eps), iters4=i; end
end

figure
hold on
plot(xx4(1,:),xx4(2,:))
bullet(xx4(1,:),xx4(2,:),15,'r')
bullet(1,1,25,'k')
contour(xcoord{1},xcoord{2},y,conts,'k:')
xlabel('x_1')
ylabel('x_2')
title('Nelder-Mead Maximization of Banana Function')
set(gca,'ytick',[0 0.5 1])
axis([-.7 1.2 -.2 1.2])

s={'Iterations';
  ['Steepest:   ' num2str(iters1)];...
  ['DFP:        ' num2str(iters2)];...
  ['BFGS:       ' num2str(iters3)];...
  ['Nelder-Mead ' num2str(iters4)]};
disp(s)

% Save Plots as EPS Files
printfigures(mfilename,4)