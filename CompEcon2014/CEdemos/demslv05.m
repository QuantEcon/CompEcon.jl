function demslv05

%% DEMSLV05 Cournot equillibrium model
%
% Solves Cournot equilibrium problem using Newton and Broyden methods

% Preliminary tasks
demosetup(mfilename)
warning off

% Model parameters
alpha = 0.6;
beta  = [0.6 0.8];

% Compute equilibrium using Newton method (explicitly)
q = [0.2;0.2];
for it=1:40
  [f,J] = cournot(q,alpha,beta);
  step = -J\f;
  q = q + step;
  if norm(step)<1.e-10, break, end
end

% Check Jacobian of cournotn
q =[1; 1];
checkjac(@cournot,q,alpha,beta);

% Generate data for contour plot
n = 100;
q1 = nodeunif(n,0.1,1.5);
q2 = nodeunif(n,0.1,1.5);
for i=1:n
  for j=1:n
    z = cournot([q1(j);q2(i)],alpha,beta);
    z1(i,j) = z(1);
    z2(i,j) = z(2);
  end
end

% Plot figure for overheads
figure
qmin = 0.1;
qmax = 1.3;
for k=1:2
  subplot(1,2,k)
  hold on
  contour(q1,q2,z1,[0 0],'k')
  contour(q1,q2,z2,[0 0],'k')
  set(gca,'xtick',[0.2 0.7 1.2])
  set(gca,'ytick',[0.2 0.7 1.2])
  q = [0.2;0.2];
  if k==1
    title('Newton''s Method')
    optset('newton','maxit',1)
    optset('newton','maxsteps',0)
    algorithm = @newton;
    func = @cournot;
  else
    title('Broyden''s Method')
    optset('broyden','maxit',1)
    optset('broyden','maxsteps',0)
    fjac = fdjac(@cournot,q,alpha,beta);
    optset('broyden','initb',inv(fjac));
    algorithm = @broyden;
    func = @cournot;
  end
  for i=1:10
    qnew = algorithm(func,q,alpha,beta);
    plot([q(1) qnew(1)],[q(2) qnew(2)],'-')
    bullet(q(1),q(2),14,'r')
    q = qnew;
  end
  axis([qmin qmax qmin qmax],'square')
  xlabel('q_1','VerticalAlignment','cap')
  ylabel('q_2','VerticalAlignment','bottom')
  ftext(0.85,qmax,'\pi''_1=0','left','top')
  ftext(qmax,0.55,'\pi''_2=0','right','middle')
end

% Save Plots as EPS Files
printfigures(mfilename,1,1,0)


%% Equilibrium function for Newton method
function [fval,fjac] = cournot(q,alpha,beta)
qsum  = q(1)+q(2);
P  = qsum^(-alpha);
P1 = (-alpha)*qsum^(-alpha-1);
P2 = (-alpha-1)*(-alpha)*qsum^(-alpha-2);
fval = [P+P1*q(1)-beta(1)*q(1); ...
        P+P1*q(2)-beta(2)*q(2)];
fjac = [2*P1+P2*q(1)-beta(1) P1+P2*q(1); ...
        P1+P2*q(2) 2*P1+P2*q(2)-beta(2)];