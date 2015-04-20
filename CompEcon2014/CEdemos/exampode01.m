% exampode01

% Script file for linear ODE examples in lecture notes.

close all
clear all


%% Basic Data
A = [-0.75 0.25; 0.25 -0.75];
x0 = [1.0;0.8];
T = 10;
  
% Layout
xlim = [-1.05 1.05];
figlayout(xlim,xlim,'Phase Diagram','x_1','x_2')
set(gca,'XTick',[-1 0 1])
set(gca,'YTick',[-1 0 1])
axis square

% Nullclines
if A(1,2)~=0&&A(2,2)~=0
  x1null = -(A(1,1)*xlim)/A(1,2);
  x2null = -(A(2,1)*xlim)/A(2,2);
  plot(xlim,x1null,'k:',xlim,x2null,'k:','LineWidth',3)
elseif A(1,2)==0
  x2null = -(A(2,1)*xlim)/A(2,2);
  plot([0 0],xlim,'k:',xlim,x2null,'k:','LineWidth',3)
elseif A(2,2)==0
  x1null = -(A(1,1)*xlim)/A(1,2);
  plot(xlim,x1null,'k:',[0 0],xlim,'k:','LineWidth',3)
end

% Eigenvalues
[V,D] = eig(A);
V = sqrt(2)*V;
disp('Eigenvalues')
disp(diag(D))
disp('Eigenvectors')
disp(V)
if isreal(D)
  for j=1:2
    if D(j,j)<0
      plot([-V(1,j) V(1,j)],[-V(2,j) V(2,j)],'g')
    else
      plot([-V(1,j) V(1,j)],[-V(2,j) V(2,j)],'b')
    end
  end
end

% Velocity Field
odevfield(@(x) A*x,xlim,xlim)

% Plot Particular Solution
t = nodeunif(200,0,T);
[V,D] = eig(A);
c = V\x0;
x = real((V*diag(c))*exp(diag(D)*t'));
bullet(x(1,1),x(2,1)+0.02,20,'r')
for j=2:200
  getframe;
  plot([x(1,j-1) x(1,j)],[x(2,j-1) x(2,j)],'r')
  if max(abs(x(:,j)))>1.05, break, end
  if max(abs(x(:,j)))<0.01, break, end
end
plot(x(1,1:j),x(2,1:j),'r')
bullet(0,0.02,20)