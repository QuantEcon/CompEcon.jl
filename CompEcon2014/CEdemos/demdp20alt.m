function demdp20alt

%% DEMDP20ALT Lifecycle Consumption-Saving Model with Pension
%
% Solves the lifecycle consumption-saving model as a nonlinear
% complementarity problem.  In this model, the agent works from period 0 to
% period T, earning a known income y(t) in each period t, a portion g of
% which he is required to pay into a pension fund. The agent retires after
% period T. At the beginning of period T+1, the agent converts the amount
% accumulated in his pension fund plus any other saving (which may be
% negative) into an annuity that provides him with fixed consumption for
% the remainder of his infinite life.  The agent begins with no net
% saving, and his saving and pension grow at a per-period rate r. The
% agent may save up to an amount slim, and may borrow up to a fixed proportion
% k of his current income.
%
% This demmo solves the Bellman equation using backward recursion.  The
% accompanying demo demdp20 solves the Euler conditions as a nonlinear
% complementarity problem using ncpsolve.
%
% Endogenous Variable
%     x       net saving by period (x<0 implies borrowing)
% Parameters
%     g       pension contribution as proportion of income
%     r       interest rate earned by assets, paid on debt
%     alpha   agent's relative risk aversion
%     slim    saving limit
%     T       number of periods until retirement
%     tmax    period of maximum income
%     ymax    maximum income (income at t=0 is 1)
%     k       borrowing limit as a proportion of income
%     delta   agent's subjective discount factor
% Derived Parameters
%     N       number of decision periods (T+1)
%     t       time index
%     R       gross interest rate (1+r)
%     y       income per period
%     Y       value of income at T+1

% Preliminary tasks
demosetup(mfilename)


%% FORMULATION
  
% Model Parameters
g     = 0.0;        % pension contribution as proportion of income
r     = 0.05;       % interest rate
alpha = 4;          % agent's relative risk aversion
slim  = 20;         % saving limit
T     = 40;         % number of periods until retirement
tmax  = 0.75*T;     % period of maximum income
ymax  = 1.5;        % maximum income (income at t=0 is 1)
k     = 1.0;        % borrowing limit as a proportion of income
delta = 0.9;        % agent's subjective discount factor

% Derived Parameters
N     = T+1;        % number of decision periods
t     = (0:T)';     % time index
R     = 1+r;        % gross interest rate
blim  = k*ymax;     % maximum borrowing limit

% Derive Income Stream
b = 2*(ymax-1)/tmax;
c = b/tmax;
y = 1 + b*t - 0.5*c*t.^2;

% Compute Current Value of Lifetime Income at T+1
Y = 0;
for i=1:N
  Y = R*(Y+y(i));
end

% Plot Income
figure
plot(t,y,'LineWidth',3)
title('Agent''s Income Stream')
xlabel('Period')
ylabel('Income')
% ylim([0 ymax])

% Define Utility Function
util = @(c) c.^(1-alpha)/(1-alpha);

% Discretize Continuous State
n    = 500;
wmin = -R*blim;
wmax =  R*slim;
w = nodeunif(n,wmin,wmax);

% Discretize Continuous Action
m = 1000;
x = nodeunif(m,-blim,slim); 

% Matrices for Discrete Optimization
W = w(:,ones(1,m));
X = x(:,ones(1,n))';


%% SOLUTION

% Initialize Value and Policy Functions
vv = zeros(n,N+1);
xx = zeros(n,N);

% Set Post-Terminal Value Function
cret = r*(w+g*Y);
U = util(cret);
U(cret<=0) = -inf;
vv(:,N+1) = U/(1-delta);

% Solve Bellman Equation
tic
for i=N:-1:1
  fprintf ('Period %5i\n',i-1)
  c = W + (1-g)*y(i) - X;
  U = util(c);
  U(c<=0|X<-k*y(i)) = -inf;
  vnext = interp1(w,vv(:,i+1),R*x);
  vnext = vnext(:,ones(1,n))';
  [vv(:,i),ix] = max(U+delta*vnext,[],2);
  xx(:,i) = x(ix);
end
toc

%% SIMULATION

% Initialize History Arrays
whist = zeros(N,1);
chist = zeros(N,1);

% Intialize
isim = 1;      % No children
wsim = 0;      % No assets

% Simulate
for i=1:N
  % Plot value function differential
  xsim = interp1(w,xx(:,i),wsim);
  csim = wsim + (1-g)*y(i) - xsim;
  whist(i) = wsim;
  chist(i) = csim;
  wsim = R*xsim;
end
% whist = [whist; wsim];

% Plot Assets & Consumption
figure
hold on
plot(t,whist,t,chist,'LineWidth',3)
legend('Assets','Consumption','Location','Best')
legend boxoff
plot(t,0*t,'k.')
title('Assets and Consumption')
xlabel('Period')

% Compute Certainty Equivalent Consumption and Income
U = sum((delta.^t).*util(chist)) + (delta^(T+1)/(1-delta))*util(r*(wsim+g*Y));
ccert = ((1-delta)*(1-alpha)*U)^(1/(1-alpha));
ycert = (r/(1+r))*sum((1/(1+r)).^t.*y);

% Print Output
fprintf('Assets at Retirement              = %5.2f\n',(wsim+g*Y))
fprintf('Retirement Income and Consumption = %5.2f\n',r*(wsim+g*Y))
fprintf('Certainty Equivalent Consumption  = %5.2f\n',ccert)
fprintf('Certainty Equivalent Income       = %5.2f\n',ycert)

%% Save Plots as EPS Files
printfigures(mfilename,2)