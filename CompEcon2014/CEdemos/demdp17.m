function demdp17

%% DEMDP17 Miscellaneous Lecture Note Figures

% Preliminary tasks
demosetup(mfilename)


%% FIGURE 1 Timber growth curve

a = 0.1;
sstar = 0.9;
c = 2*(sstar-a)/sstar^2;
b = c*sstar;

% Define timber growth curve
F = @(s) a+b*s-0.5*c*s.^2;

% Figure setup
n = 401;
smin = 0;
smax = 1;
s  = nodeunif(n,smin,smax);

% Plot timber growth curve
figure
hold on
% Figure setup
axis square
xlim ([0 1])
ylim ([0 1])
set(gca,'xtick',[])
set(gca,'ytick',[])
% Plot growth function
plot(s,F(s))
% Plot 45-degree line
plot(s,s,'k--','Linewidth',2)
ftext(0.12,0.01,'45^\circ','right','bottom',10);
% Plot fixed points
plot(sstar,sstar,'k*','linewidth',8)
plot([sstar sstar],[0 sstar],'k--','linewidth',2)
ftext(0,0.95,'h(s)','right','middle',14);
ftext(1,0-.015,'s','center','top',14);
ftext(0.9,0-.015,'s','center','top',14);
ftext(0.9-.006,0.01,'\bf -','center','top',14);
ftext(0-0.01,a,'h(0)','right','middle',12);


%% FIGURE 2 Asset productivity curves

% Figure setup
n = 401;
amin = 0;
amax = 1;
a  = nodeunif(n,amin,amax);

% Define timber growth curve
f1 = @(a) 0.9-0.4*(exp(a)-1);
f2 = @(a) 0.9-2.7*(a-0.5).^2;
f1(0)

figure
subplot(1,2,1)
hold on
plot(a,f1(a))
xlabel('Age','FontSize',16)
ylabel('Output','FontSize',16)
title('Declining','FontSize',16)
set(gca,'xtick',[])
set(gca,'ytick',[])
xlim ([0 1])
ylim ([0 1])
axis square
subplot(1,2,2)
hold on
plot(a,f2(a))
xlabel('Age','FontSize',16)
ylabel('Output','FontSize',16)
title('Concave','FontSize',16)
set(gca,'xtick',[])
set(gca,'ytick',[])
xlim ([0 1])
ylim ([0 1])
axis square


%% FIGURE 3 Consumer Surplus

% Initiate figure
figure
hold on

% Define marginal cost curve
f = @(q) 0.2+0.7*exp(-5*q);

% Set plotting parameters
n = 1001;
qmin = 0;
qmax = 1;
pmin = 0;
pmax = 1;
q = nodeunif(n,0,qmax);
p = f(q);

% Plot area under inverse demand curve
qstar = 0.4;
pstar = f(qstar);
kstar = 0.5*pstar;
for i=1:n
  if 0.005<q(i) && q(i)<qstar
    plot([q(i) q(i)],[kstar p(i)-0.02],'y-','linewidth',2)
  end
end

% Plot inverse demand curve
plot(q,p)
xlim([qmin qmax])
ylim([pmin pmax])
set(gca,'xtick',[])
set(gca,'ytick',[])
axis square
% xlabel('Quantity')
% ylabel('Price')
ftext(qmin,pmin,'0','right','top');
ftext(qmax,pmin,'q','center','top');
ftext(qmin-0.02,pmax,'p','right','top');

% Annotate figure
ftext(qmin-0.02,pstar,'p_t','right','middle');
ftext(qmin-0.02,kstar,'k','right','middle');
ftext(qstar,pmin,'q_t','center','top');
plot([0 qstar],[pstar pstar],'k--','linewidth',2)
plot([0 qstar],[kstar kstar],'k--','linewidth',2)
plot([qstar qstar],[pmin pstar],'k--','linewidth',2)
qtmp = 0.4*qstar; ptmp = 0.5*f(qtmp);
ftext(qtmp,ptmp-0.04,'PS','center','middle');
ftext(qtmp,ptmp+0.11,'CS','center','middle');


%% FIGURE 4 Renewable resource model biological growth function

% Gross return
R  = 1.2;

% Define biological growth function
m1 = 0.4;
m2 = 0.6;
s1 = 1.0;
s2 = 1.0;
c  = 0.45;
F1 = @(q) cdf('norm',q,m1,(s1*m1)^2);
F2 = @(q) cdf('norm',q,m2,(s2*m2)^2);
f1 = @(q) pdf('norm',q,m1,(s1*m1)^2);
f2 = @(q) pdf('norm',q,m2,(s2*m2)^2);
c2 = c; c1 = c*(f2(0))/(f1(0));
F  = @(q) c1*(F1(q)-F1(0))-c2*(F2(q)-F2(0));

% Compute fixed points
G  = @(q) F(q)-q;
q1 = broyden(G,m1);
q2 = broyden(G,m2);

% Compute optimal retention, stock, harvest
g  = @(q) c1*f1(q)-c2*f2(q)-R;
rstar = broyden(g,m2);
sstar = F(rstar);

% Figure setup
n = 401;
qmin = 0;
qmax = 1;
q = nodeunif(n,qmin,qmax);

% Plot growth function and steady-states
figure
hold on
% Figure setup
axis square
ylim ([0 1])
set(gca,'xtick',[])
set(gca,'ytick',[])
% Plot growth function
plot(q,F(q))
% Plot 45-degree line
plot(q,q,'k--','Linewidth',2)
ftext(0.08,0.08,'45^\circ','right','bottom',10);
% Plot fixed points
plot(0 ,0 ,'k*','linewidth',8)
plot(q1,q1,'k*','linewidth',8)
plot(q2,q2,'k*','linewidth',8)
plot([q1 q1],[0 q1],'k--','linewidth',2)
plot([q2 q2],[0 q2],'k--','linewidth',2)
ftext(q1+0.05,0.1,'s^*_1','center','top',14);
ftext(q2+0.05,0.1,'s^*_2','center','top',14);
xlabel('Stock This Period')
ylabel('Stock Next Period')


%% FIGURE 5 Renewable resource model optimal solution

% Plot growth function and optimal solution
figure
hold on
% Figure setup
axis square
ylim ([0 1])
set(gca,'xtick',[])
set(gca,'ytick',[])
% Plot growth function
plot(q,F(q))
% Plot 45-degree line
plot(q,q,'k--','Linewidth',2)
ftext(0.08,0.08,'45^\circ','right','bottom',10);
% Plot slope line
qq = nodeunif(n,rstar-0.2,rstar+0.2);
yy = F(rstar)+R*(qq-rstar);
plot(qq,yy,'k--','linewidth',2)
plot(rstar,sstar,'k*','linewidth',8)
% Plot optimal retention, stock, harevst
fudge = 0.02;
plot([0 rstar],[sstar sstar],'k--','linewidth',2)
plot([0 rstar],[rstar rstar],'k--','linewidth',2)
plot([rstar rstar],[0 rstar],'k--','linewidth',2)
plot([rstar rstar],[rstar sstar],'k--','linewidth',2)
ftext(0-fudge,rstar,'s^*-q^*','right','middle',14);
ftext(0-fudge,sstar,'s^*','right','middle',14);
plot(rstar,rstar,'k*','linewidth',8)
ftext(0.65,0.9,'slope=1+\rho','right','bottom',10);
ftext(rstar,0,'s^*-q^*','center','top',14);
mid = 0.5*(rstar+sstar);
ftext(0+fudge,mid,'q^*','left','middle',14);
annotation('arrow',[0.23 0.23],[mid-0.04 rstar+0.01],'LineWidth',2)
annotation('arrow',[0.23 0.23],[mid+0.01 sstar-0.04],'LineWidth',2)


%% FIGURE 6 Cost of Extraction 

% Initiate figure
figure
hold on

% Define marginal cost curve
f = @(s) 0.2+0.7*exp(-5*s);

% Set plotting parameters
n = 401;
smin = 0;
smax = 1;
kmin = 0;
kmax = 1;
s = nodeunif(n,smin,smax);

% Plot marginal cost curve
k = f(s);
plot(s,k)
xlim([smin smax])
ylim([kmin kmax])
set(gca,'xtick',[])
set(gca,'ytick',[])
% xlabel('Ore Stock')
% ylabel('Marginal Cost of Extraction')
ftext(smin,kmin,'0','right','top');
ftext(smin,kmax,'k(s)','right','top');
ftext(smax,f(smax)-0.2,'s','center','top');

% Plot area under marginal cost curve
a = 0.50;
b = 0.75;
for i=1:n
  if a<=s(i) && s(i)<=b
    plot([s(i) s(i)],[kmin k(i)-0.01],'y-','linewidth',2)
  end
end
plot([a a],[kmin f(a)],'k--','linewidth',2)
plot([b b],[kmin f(b)],'k--','linewidth',2)
ftext(a,kmin,'s_t-q_t','center','top');
ftext(b,kmin,'s_t','center','top');
stmp = 0.5*(a+b); ktmp = 0.5*f(stmp);
ftext(stmp,ktmp,'C_t','center','middle');

% Pot abandonment point
sstar = 0.25;
kstar = f(sstar);
ftext(sstar,kmin,'s^*','center','top');
plot([sstar sstar],[kmin kstar],'k--','linewidth',2)
plot([smin sstar],[kstar kstar],'k--','linewidth',2)
ftext(smin,kstar,'p(0)','right','middle');


%% Save Plots as EPS Files
printfigures(mfilename,6,1,0)