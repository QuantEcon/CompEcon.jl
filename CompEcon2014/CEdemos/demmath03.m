function demmath03

%% DEMMATH03 Major Continuous Distribution CDFs and PDFs

% Preliminary tasks
demosetup(mfilename)


%% Normal Distribution

n = 400;
xmin = -4;
xmax =  4;
x = nodeunif(n,xmin,xmax);

mu  = 0;
var = 1;

layout([-4;0;4],[],'pdf')
y = pdf('Normal',x,mu,var);
plot(x,y,'b','LineWidth',4)

layout([-4;0;4],[0;1],'cdf')
y = cdf('Normal',x,mu,var);
xcrit = icdf('Normal',0.95,mu,var)
plot([xcrit xcrit],[0 0.95],'r:','LineWidth',4)
plot([-4 xcrit],[0.95 0.95],'r:','LineWidth',4)
text(-4,0.95,'\bf 0.95','VerticalAlignment','middle','HorizontalAlignment','right','Fontsize',12)
text(xcrit,0,'\bf x_{0.95}','VerticalAlignment','top','HorizontalAlignment','center','Fontsize',12)
plot(x,y,'b','LineWidth',4)
ylim([0 1])

prob = cdf('Normal',1,mu,var)-cdf('Normal',-1,mu,var)
prob = cdf('Normal',2,mu,var)-cdf('Normal',-2,mu,var)
xcrit = icdf('Normal',0.950,mu,var)
xcrit = icdf('Normal',0.975,mu,var)


%% Lognormal Distribution

n = 400;
x = nodeunif(n,0,5);

mu  = 0;
var = 1;

layout([0;5],[],'pdf')
y = pdf('Lognormal',x,mu,var);
plot(x,y,'b','LineWidth',4)

layout([0;5],[0;1],'cdf')
y = cdf('Lognormal',x,mu,var);
plot(x,y,'b','LineWidth',4)
xcrit = icdf('Lognormal',0.90,mu,var)
plot([xcrit xcrit],[0 0.9],'r:','LineWidth',4)
plot([0 xcrit],[0.9 0.9],'r:','LineWidth',4)
ylim([0 1])
text(0,0.9,'\bf 0.9','VerticalAlignment','middle','HorizontalAlignment','right','Fontsize',12)
text(xcrit,0,'\bf x_{0.9}','VerticalAlignment','top','HorizontalAlignment','center','Fontsize',12)

mu  = 1;
var = 2;
x90 = icdf('Lognormal',0.90,mu,var)
x95 = icdf('Lognormal',0.95,mu,var)
x99 = icdf('Lognormal',0.99,mu,var)


%% Beta Distribution

n = 400;
x = nodeunif(n,0,5);

beta = 1;

layout([0;5],[],'pdf')
y = pdf('Exponential',x,beta);
plot(x,y,'b','LineWidth',4)

layout([0;5],[0;1],'cdf')
y = cdf('Exponential',x,beta);
xmode = icdf('Exponential',0.5,beta)
plot(x,y,'b','LineWidth',4)
plot([xmode xmode],[0 0.5],'r:','LineWidth',4)
plot([0 xmode],[0.5 0.5],'r:','LineWidth',4)
ylim([0 1])
text(0,0.5,'\bf 0.5','VerticalAlignment','top','HorizontalAlignment','center','Fontsize',12)
text(xmode,0,'\bf x_{0.5}','VerticalAlignment','top','HorizontalAlignment','center','Fontsize',12)


%% Gamma Distribution

n = 400;
x = nodeunif(n,0,5);

layout([0;5],[],'pdf')
alpha = 1;
theta  = 0.5;
y = pdf('Gamma',x,alpha,theta);
plot(x,y,'b','LineWidth',4)
alpha = 2;
theta  = 0.5;
y = pdf('Gamma',x,alpha,theta);
plot(x,y,'r','LineWidth',4)
legend('\alpha=1,\theta=0.5','\alpha=2,\theta=0.5','Location','ne')
legend boxoff

layout([0;5],[0;1],'cdf')
alpha = 1;
theta  = 0.5;
y = cdf('Gamma',x,alpha,theta);
plot(x,y,'b','LineWidth',4)
alpha = 2;
theta  = 0.5;
y = cdf('Gamma',x,alpha,1/beta);
plot(x,y,'r','LineWidth',4)
ylim([0 1])
legend('\alpha=1,\theta=0.5','\alpha=2,\theta=0.5','Location','e')
legend boxoff


%% Beta Distribution

n = 400;
x = nodeunif(n,0,1);

layout([0;1],[],'pdf')
alpha = 3;
beta  = 3;
y = pdf('Beta',x,alpha,beta);
plot(x,y,'b','LineWidth',4)
alpha = 2;
beta  = 5;
y = pdf('Beta',x,alpha,beta);
plot(x,y,'r','LineWidth',4)
alpha = 5;
beta  = 2;
y = pdf('Beta',x,alpha,beta);
plot(x,y,'g','LineWidth',4)
legend('\alpha=3,\beta=3','\alpha=2,\beta=5','\alpha=5,\beta=2','Location','n')
ylim([0 3])
legend boxoff

layout([0;1],[0;1],'cdf')
alpha = 3;
beta  = 3;
y = cdf('Beta',x,alpha,beta);
plot(x,y,'b','LineWidth',4)
alpha = 2;
beta  = 5;
y = cdf('Beta',x,alpha,beta);
plot(x,y,'r','LineWidth',4)
alpha = 5;
beta  = 2;
y = cdf('Beta',x,alpha,beta);
plot(x,y,'g','LineWidth',4)
ylim([0 1])
legend('\alpha=3,\beta=3','\alpha=2,\beta=5','\alpha=5,\beta=2','Location','nw')
legend boxoff


%% Chi-Squared Distribution

n = 400;
x = nodeunif(n,0,15);

layout([0;12],[],'pdf')
k = 2;
y = pdf('Chisquare',x,k);
plot(x,y,'b','LineWidth',4)
k = 3;
y = pdf('Chisquare',x,k);
plot(x,y,'r','LineWidth',4)
k = 5;
y = pdf('Chisquare',x,k);
plot(x,y,'g','LineWidth',4)
legend('k=2','k=3','k=5','Location','n')
ylim([0 0.45])
legend boxoff

layout([0;12],[0;1],'cdf')
k = 2;
y = cdf('Chisquare',x,k);
plot(x,y,'b','LineWidth',4)
k = 3;
y = cdf('Chisquare',x,k);
plot(x,y,'r','LineWidth',4)
k = 5;
y = cdf('Chisquare',x,k);
plot(x,y,'g','LineWidth',4)
legend('k=2','k=3','k=5','Location','e')
ylim([0 1])
legend boxoff


%% F Distribution

n = 400;
x = nodeunif(n,0,5);

layout([0;5],[],'pdf')
m = 1;
n = 2;
y = pdf('F',x,n,m);
plot(x,y,'b','LineWidth',4)
m = 2;
n = 5;
y = pdf('F',x,n,m);
plot(x,y,'r','LineWidth',4)
m = 100;
n = 10;
y = pdf('F',x,n,m);
plot(x,y,'g','LineWidth',4)
legend('n=2,m=1','n=5,m=2','n=100,m=10','Location','n')
legend boxoff

layout([0;5],[0;1],'cdf')
m = 1;
n = 2;
y = cdf('F',x,n,m);
plot(x,y,'b','LineWidth',4)
m = 2;
n = 5;
y = cdf('F',x,n,m);
plot(x,y,'r','LineWidth',4)
m = 100;
n = 10;
y = cdf('F',x,n,m);
plot(x,y,'g','LineWidth',4)
ylim([0 1])
legend('n=2,m=1','n=5,m=2','n=100,m=10','Location','se')
legend boxoff


%% Student's T Distribution

n = 400;
x = nodeunif(n,-4,4);

layout([-4;0;4],[],'pdf')
k = 1;
y = pdf('T',x,k);
plot(x,y,'b','LineWidth',4)
k = 4;
y = pdf('T',x,k);
plot(x,y,'r','LineWidth',4)
k = 100;
y = pdf('T',x,k);
plot(x,y,'g','LineWidth',4)
legend('\nu=1','\nu=4','\nu=100','Location','ne')
legend boxoff

layout([-4;0;4],[0;1],'cdf')
k = 1;
y = cdf('T',x,k);
plot(x,y,'b','LineWidth',4)
k = 4;
y = cdf('T',x,k);
plot(x,y,'r','LineWidth',4)
k = 100;
y = cdf('T',x,k);
plot(x,y,'g','LineWidth',4)
legend('\nu=1','\nu=4','\nu=100','Location','se')
legend boxoff


%% Logistic Distribution

n = 400;
xmin = -8;
xmax =  8;
x = nodeunif(n,xmin,xmax);

mu  = 0;
sigma = 1;

layout([xmin;0;xmax],[],'pdf')
y = pdf('Logistic',x,mu,sigma);
plot(x,y,'b','LineWidth',4)
xlim([xmin xmax])

layout([xmin;0;xmax],[0;1],'cdf')
y = cdf('Logistic',x,mu,sigma);
xcrit = icdf('Logistic',0.9,mu,sigma)
plot([xcrit xcrit],[0 0.9],'r:','LineWidth',4)
plot([xmin xcrit],[0.9 0.9],'r:','LineWidth',4)
text(xmin,0.9,'\bf 0.9','VerticalAlignment','middle','HorizontalAlignment','right','Fontsize',12)
text(xcrit,0,'\bf x_{0.9}','VerticalAlignment','top','HorizontalAlignment','center','Fontsize',12)
plot(x,y,'b','LineWidth',4)
xlim([xmin xmax])
ylim([0 1])

prob  = cdf('Logistic',1,mu,sigma)-cdf('Logistic',-1,mu,sigma)
prob  = cdf('Logistic',2,mu,sigma)-cdf('Logistic',-2,mu,sigma)
xcrit = icdf('Logistic',0.950,mu,sigma)
xcrit = icdf('Logistic',0.975,mu,sigma)


% Additional Plots

figure
n = 3;
p = 0.6;
x = 0:n; 
y = pdf('Binomial',x,n,p);
bar(x,y,0.8)
set(gca,'XTick',0:3,'Fontsize',14)
set(gca,'YTick',0:0.1:0.5,'Fontsize',14)
xlabel('x','Fontsize',14)
ylabel('Probability Mass','Fontsize',14)

figure
n = 6;
p = 0.6;
x = 0:n; 
y = pdf('Geometric',x,p);
bar(x,y,0.8)
xlim([-0.5,n+0.5])
set(gca,'XTick',0:n,'Fontsize',14)
set(gca,'YTick',0:0.1:0.6,'Fontsize',14)
xlabel('x','Fontsize',14)
ylabel('Probability Mass','Fontsize',14)

figure
n = 8;
lambda = 1.5;
x = 0:n; 
y = pdf('Poisson',x,lambda);
bar(x,y,0.8)
xlim([-0.5,n+0.5])
set(gca,'XTick',0:n,'Fontsize',14)
set(gca,'YTick',0:0.1:0.3,'Fontsize',14)
xlabel('x','Fontsize',14)
ylabel('Probability Mass','Fontsize',14)





%% Save Plots as EPS Files
printfigures(mfilename,21)


%% Layout Function
function layout(xticks,yticks,title)
figure
hold on
box off
xlabel('x','Fontsize',14)
if strcmp(title,'pdf')
    ylabel('Probability Density','Fontsize',14)
else
    ylabel('Cumulative Probability','Fontsize',14)
end
set(gca,'XTick',xticks,'Fontsize',14)
set(gca,'YTick',yticks,'Fontsize',14)