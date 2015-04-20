function demode01

%% DEMODE01 - Stability of Linear Homogeneous ODEs

% Preliminary tasks
demosetup(mfilename)


%% Example 1
A = [-0.75 0.25; 0.25 -0.75];
x0 = [1.0;0.8];
T = 10;
phase(A,'Negative Real Eigenvalues')
plotnullcline(-0.65,-0.98,-1.00,-0.45)
ploteigenvector(-0.95,1.00,0.82,1.00)
path(200,T,x0,A);


%% Example 2
A = [0.75 -0.25; -0.25 0.75];
x0 = [-0.1;-0.2];
T = 3;
phase(A,'Positive Real Eigenvalues')
plotnullcline(-0.65,-0.98,-1.00,-0.45)
ploteigenvector(0.82,1.00,-0.95,1.00)
path(200,T,x0,A);


%% Example 3
A = [-0.25 0.75; 0.75 -0.25];
x0 = [1.0;-0.8];
T = 5;
phase(A,'Real Eigenvalues of Oposite Signs')
plotnullcline(-1.00,-0.45,-0.65,-0.98)
ploteigenvector(-0.95,1.00,0.82,1.00)
path(200,T,x0,A);


%% Example 4
A = [-0.5 0; 0 -0.5];
x0 = [0.5;0.5];
T = 8;
phase(A,'One Repeated Negative Real Eigenvalue')
plotnullcline(-0.32,-0.98,-1.00,-0.12)
path(200,T,x0,A);


%% Example 5
A = [0.5 0; 0 0.5];
x0 = [0.5;0.5];
T = 2;
phase(A,'One Repeated Positive Real Eigenvalue')
plotnullcline(-0.32,-0.98,-1.00,-0.12)
path(200,T,x0,A);


%% Example 6
A = [0 1;-1 0];
x0 = [0.5;0.5];
T = 8;
phase(A,'Complex Eigenvalues with Zero Real Part')
plotnullcline(-1.00,-0.12,-0.32,-0.98)
path(200,T,x0,A);


%% Example 7
A = [0 1;-1 -0.3];
x0 = [0.7;0.7];
T = 30;
phase(A,'Complex Eigenvalues with Negative Real Part')
plotnullcline(-1.00,-0.12,0.00,-0.98)
path(200,T,x0,A);


%% Example 8
A = [0 1;-1 0.1];
x0 = [0.2;0.2];
T = 30;
phase(A,'Complex Eigenvalues with Positive Real Part')
plotnullcline(-1.00,-0.12,-0.42,-0.98)
path(200,T,x0,A);


%% Save Plots as EPS Files
printfigures(mfilename,8)


function phase(A,figtitle)

% Layout
xlim = [-1.05 1.05];
figlayout(xlim,xlim,figtitle,'x_1','x_2')
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
odefield(@(x) A*x,xlim,xlim)


%% Plot Particular Solution
function path(N,T,x0,A)
t = nodeunif(N,0,T);
[V,D] = eig(A);
c = V\x0;
x = real((V*diag(c))*exp(diag(D)*t'));
bullet(x(1,1),x(2,1)+0.02,20,'r')
for j=2:N
  getframe;
  plot([x(1,j-1) x(1,j)],[x(2,j-1) x(2,j)],'r')
  if max(abs(x(:,j)))>1.05, break, end
  if max(abs(x(:,j)))<0.01, break, end
end
plot(x(1,1:j),x(2,1:j),'r')
bullet(0,0.02,20)
% disp('Hit Return to Continue to Next Example')
% pause


%% Plot Nullcline Labels
function plotnullcline(x11,x12,x21,x22)
text(x11,x12,'x''_1=0','FontSize',12)
text(x21,x22,'x''_2=0','FontSize',12)


%% Plot Eigenvalue Labels
function ploteigenvector(x11,x12,x21,x22)
text(x11,x12,'v_1','FontSize',12)
text(x21,x22,'v_2','FontSize',12)