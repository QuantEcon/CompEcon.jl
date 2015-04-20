function demapp07

%% DEMAPP07 Solve Cournot oligopoly model via collocation

% Preliminary tasks
demosetup(mfilename)

% Model parameters
alpha = 1.0;
eta   = 3.5;


%% Residual function
resid = @(c,basis,p,alpha,eta) ...
  p + funeval(c,basis,p).*((-1./eta)*p.^(eta+1)) ...
    - alpha*sqrt(funeval(c,basis,p)) ...
    - funeval(c,basis,p).^2;

% Approximation structure
n =  21;
a = 0.5;
b = 2.0;
basis = fundefn('cheb',n,a,b);
pnode = funnode(basis);

% Solve for effective supply function
c = [1;zeros(n-1,1)];
optset('broyden','tol',10e-12);
c = broyden(resid,c,basis,pnode,alpha,eta);

% Plot demand and effective supply for m=5 firms
figure
pplot = nodeunif(501,a,b);
splot = funeval(c,basis,pplot);
dplot = pplot.^-eta;
plot(5*splot,pplot,dplot,pplot);
xlim([0 4])
ylim([0.5 2])
legend('Supply','Demand')
title('Cournot Effective Firm Supply Function')
xlabel('Quantity'); 
ylabel('Price');

% Plot residual
figure
hold on
rplot = resid(c,basis,pplot,alpha,eta);
plot(pplot,rplot)
plot(pplot,zeros(size(rplot)),'k--','LineWidth',2)
title('Residual Function for Cournot Problem')
xlabel('Price'); 
ylabel('Residual')

% Plot demand and effective supply for m=1,3,5,10,15,20 firms
figure
m = [1 3 5 10 15 20];
plot(splot*m,pplot,pplot.^(-eta),pplot)
title('Industry Supply and Demand Functions')
xlabel('Quantity'); 
ylabel('Price')
legend('m=1','m=3','m=5','m=10','m=15','m=20');
xlim([0 13]);

% Plot equilibrium price as a function of number of firms m
figure
pp = (b+a)/2;
dp = (b-a)/2;
m  = (1:25)';
for i=1:50
    dp = dp/2;
    pp = pp-sign(funeval(c,basis,pp).*m-pp.^(-eta)).*dp;
end
plot(m,pp)
title('Cournot Equilibrium Price as Function of Industry Size')
xlabel('Number of Firms'); 
ylabel('Price')

% Save Plots as EPS Files
printfigures(mfilename,4)