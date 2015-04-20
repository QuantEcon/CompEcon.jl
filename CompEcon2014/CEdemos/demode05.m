%% DEMODE05 Commodity Storage Model
%
%  Solve
%    s' = -p^(-eta)
%    p' = r*p+k
%  where
%    s: stocks
%    p: price
%    x = [s;p]

% Preliminary tasks
demosetup(mfilename)


%% FORMULATION

% Model Parameters
r   = 0.1;                              % interest rate
k   = 0.5;                              % unit cost of storage
eta = 5;                                % demand elasticity

% Velocity Function
f = @(x) [-x(2,:).^(-eta);r*x(2,:)+k];


%% SOLVE ODE USING ODECOL

% Solve ODE
n   = 15;                     % number of basis functions
T   = 1;                      % time horizon
s0  = 1;                      % initial stocks
sT  = 0;                      % terminal stocks
bx = [1;1];                   % boundary variables
bt = [0;T];                   % boundary times
bv = [s0;sT];                 % boundary values
c = zeros(n,2); c(1,:) = 1;
[x,t,res] = odecol(f,bv,T,n,bt,bx,[],c);


%% PLOT SOLUTION

% Plot Solution
figure
plot(t,x)
title('Commodity Storage Model')
ylim([0 1.5])
xlabel('Time')
legend('Stocks','Price')
set(legend,'Box','off')
set(legend,'Location','NorthWest')
set(legend,'FontSize',14)

% Plot Residuals
figure
plot(t,res,t,0*res,'k:')
title('Commodity Storage Model')
xlabel('Time')
ylabel('Residual')
legend('Stocks','Price')
set(legend,'Box','off')
set(legend,'Location','NorthWest')
set(legend,'FontSize',14)


%% Save Plots as EPS Files
printfigures(mfilename,2)