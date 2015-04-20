function demdp20

%% DEMDP20 Lifecycle Consumption-Saving Model with Pension
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
% agent may save up to an amount L, and may borrow up to a fixed proportion
% k of his current income.
%
% This demo solves the Euler conditions as a nonlinear complementarity
% problem using ncpsolve.  The accompanying demo demdp20alt solves the
% Bellman equation using backward recursion.
%
% Endogenous Variable
%     x       net saving by period (x<0 implies borrowing)
% Parameters
%     g       pension contribution as proportion of income
%     r       interest rate earned by assets, paid on debt
%     alpha   agent's relative risk aversion
%     L       saving limit
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
r     = 0.05;       % interest rate earned by assets, paid on debt
alpha = 4;          % agent's relative risk aversion
L     = 20;         % saving limit
T     = 40;         % number of periods until retirement
tmax  = 0.75*T;     % period of maximum income
ymax  = 1.5;        % maximum income (income at t=0 is 1)
k     = 1.0;        % borrowing limit as a proportion of income
delta = 0.9;        % agent's subjective discount factor

% Derived Parameters
N     = T+1;        % number of decision periods
t     = (0:T)';     % time index
R     = 1+r;        % gross interest rate

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

% Define Utility and Derivatives
u  = @(c) c.^(1-alpha)/(1-alpha);
u1 = @(c) c.^(-alpha);
u2 = @(c) -alpha*c.^(-alpha-1);


%% SOLUTION

% Solve Euler Conditions
x = 0.1*y;
[x,e] = ncpsolve(@F,-k*y,L,x,y,delta,R,r,g,Y,u1,u2);

% Check Function File for Internal Consistency
% [error,i,j] = checkjac(@F,x,y,delta,R,r,g,Y,u1,u2);

% Compute Consumption
c = R*[0;x(1:N-1,1)] + (1-g)*y - x;
if any(c<0), disp('WARNING: Negative Consumption'), end

% Compute Assets
s = R*[0;x(1:N-1,1)];
if any(c<0), disp('WARNING: Negative Consumption'), end

% Plot Assets & Consumption
figure
hold on
plot(t,s,t,c,'LineWidth',3)
plot(t,-k*y,'r:','LineWidth',2)
legend('Assets','Consumption','Borrowing Limit','Location','Best')
legend boxoff
plot(t,0*t,'k.')
title('Assets and Consumption')
xlabel('Period')

% Compute Certainty Equivalent Consumption and Income
U = sum((delta.^t).*u(c)) + (delta^(T+1)/(1-delta))*u(R*r*(g*Y+x(N)));
ccert = ((1-delta)*(1-alpha)*U)^(1/(1-alpha));
ycert = (r/(1+r))*sum((1/(1+r)).^t.*y);

% Print Output
fprintf('Assets at Retirement              = %5.2f\n',R*(g*Y+x(N)))
fprintf('Retirement Income and Consumption = %5.2f\n',R*r*(g*Y+x(N)))
fprintf('Certainty Equivalent Consumption  = %5.2f\n',ccert)
fprintf('Certainty Equivalent Income       = %5.2f\n',ycert)
fprintf('Norm of Euler Equation            = %9.2e\n',norm(e))


%% PARAMETRIC ANALYSIS

% Solve for Different Borrowing Limits
kvec = [0 1 2];

% Solve Euler Conditions
xinit = 0.1*y;
x = zeros(N,length(kvec));
for i=1:length(kvec)
  x(:,i) = ncpsolve(@F,-kvec(i)*y,L,xinit,y,delta,R,r,g,Y,u1,u2);
end
s = R*[zeros(1,length(kvec));x(1:N-1,:)];

% Plot Assets
figure
hold on
plot(t,s,'LineWidth',3)
legend('k=0','k=1','k=2','Location','Best')
legend boxoff
plot(t,0*t,'k.')
title('Assets - Different Borrowing Limits')
xlabel('Period')
ylabel('Assets')

% Solve for Different Interest Rates
rvec = [0.05 0.10 0.15];
Rvec = 1+rvec;

% Solve Euler Conditions
xinit = 0.1*y;
x = zeros(N,length(rvec));
for i=1:length(rvec)
  Ytmp = 0;
  for j=2:N
    Ytmp = Rvec(i)*Ytmp+y(j);
  end
  x(:,i) = ncpsolve(@F,-k*y,L,xinit,y,delta,Rvec(i),rvec(i),g,Ytmp,u1,u2);
end
s = R*[zeros(1,length(rvec));x(1:N-1,:)];

% Plot Assets
figure
hold on
plot(t,s,'LineWidth',3)
legend('r=0.05','r=0.10','r=0.15','Location','Best')
legend boxoff
plot(t,0*t,'k.')
title('Assets - Different Interest Rates')
xlabel('Period')
ylabel('Assets')

% Solve for Different Pension Contributions
gvec = [0.0 0.1 0.2];

% Solve Euler Conditions
xinit = 0.1*y;
x = zeros(N,length(rvec));
for i=1:length(rvec)
  x(:,i) = ncpsolve(@F,-k*y,L,xinit,y,delta,R,r,gvec(i),Y,u1,u2);
end
s = R*[zeros(1,length(gvec));x(1:N-1,:)];

% Plot Assets
figure
hold on
plot(t,s,'LineWidth',3)
legend('g=0.0','g=0.1','g=0.2','Location','Best')
legend boxoff
plot(t,0*t,'k.')
title('Assets - Different Pension Fund Contributions')
xlabel('Period')
ylabel('Assets')


%% Save Plots as EPS Files
printfigures(mfilename,5)


%% Function File
%
function [f,J] = F(x,y,delta,R,r,g,Y,u1,u2)
N = length(x);
c = R*[0;x(1:N-1,1)] + (1-g)*y - x;
f = -u1(c) + [delta*R*  u1(c(2:N,1)); (delta/(1-delta))*R*  r*  u1(R*r*(g*Y+x(N)))];
J =  u2(c) + [delta*R^2*u2(c(2:N,1)); (delta/(1-delta))*R^2*r^2*u2(R*r*(g*Y+x(N)))];
J = diag(J);
for i=1:N-1
  J(i,i+1) = -delta*R*u2(c(i+1));
end
for i=2:N
  J(i,i-1) = -R*u2(c(i));
end
%%