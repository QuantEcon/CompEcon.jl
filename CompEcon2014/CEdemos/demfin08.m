function demfin08

%% DEMFIN08 Affine Asset Pricing

% Preliminary tasks
demosetup(mfilename)


%% FORMULATION

% Model Parameters
T      =  30;         % time-to-maturity
kappa  =  0.1;        % speed of mean reversion
alpha  =  0.05;       % long-run mean interest rate 
sigma  =  0.1;        % interest rate volatility

% Convert parameters to standard form
a = kappa*alpha;
A = -kappa;
C = sigma;
b = 0;
B = 1;


%% SOLUTION

% Solve Model
tau = linspace(0,T,301)';
[beta,beta0] = affasset(tau,a,A,b,B,C,1,0,0,0);


%% ANALYSIS

% Create plots
r = linspace(0,0.25,101)';
V = exp(beta0(end)+beta(end)*r);

figure
plot(r,V)
title([num2str(T) ' Year Zero-Coupon Bond Price'])
xlabel('r')
ylabel('V(r)')

figure
plot(tau,[beta beta0])
title('\beta and \beta_0')
xlabel('Time to Maturity')

r = (0.03:0.01:0.08);
m = length(r);

R = -(beta0(:,ones(m,1))+beta*r)./tau(:,ones(m,1));
R(1,:) = r;
figure
plot(tau,R)
xlabel('Time to Maturity')
ylabel('Yield')
title('Term Structures for Alternative Short Rates')

nn = length(r);
legend([repmat('r  =  ',nn,1) reshape(sprintf('%4.2f',r),4,nn)'])

% Save Plots as EPS Files
printfigures(mfilename,3)