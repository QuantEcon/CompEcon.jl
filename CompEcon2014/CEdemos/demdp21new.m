function demdp21new

%% WARNING: UNDER DEVELOPMENT WITH TONY GALLENSTEIN.


%% Sterilization Model

% States
%   w   wealth (w<0 implies indebted)
%   i   number of children plus 1
% Action 
%   x   savings (x<0 implies borowing)
%   j   decision to have child (1=no, 2=yes)

% Preliminary tasks
close all
clear all


%% FORMULATION

% Model Parameters
T1    = 15;             % fertility decision periods
T2    = 40;             % savings decision periods
Y     = 1;              % permanent income
blim  = 0.40;           % borrowing limit
slim  = 6.00;           % saving limit
k     = 0.00;           % cost of contraception
p	  = 0.05;           % net cost per child
P     = 2.50;           % payoff per child at retirement
g     =	zeros(T1+1,1);  % monetized benefit from children
cmin  = 0.05;           % minimum consumption
alpha = 2.0;            % relative risk aversion
r     = 0.05;           % interest rate
R     = 1+r;            % gross interest rate
delta = 0.90;           % subjective discount rate
sigma = 0.002;            % transitory income shock volatility

% Transitory Income Shock Distribution
nshk = 5;                                 	% number of income shocks
[e,q] = qnwlogn(nshk,-sigma^2/2,sigma^2); 	% income shocks and probabilities
y = Y*e;                                    % incomes

% Utility Function
util = @(c) c.^(1-alpha)/(1-alpha);

% Discretize Wealth
n    = 1000;
wmin = -R*blim+min(y)
wmax =  R*slim+max(y)
w = nodeunif(n,wmin,wmax);

% Discretize Net Saving
m = 500;
x = nodeunif(m,-blim,slim); 

% Matrices for Fast Discrete Optimization
W    = w(:,ones(1,m));
X    = x(:,ones(1,n))';
Wnxt = R*x(:,ones(1,nshk))+y(:,ones(1,m))';

% Compute Utilities
U0 = zeros(n,m,T1+1);   % do not use contraception
U1 = zeros(n,m,T1+1);   % use contraception
for i=1:T1+1
  C = W - X - p*(i-1) + g(i);
  Utmp = util(C);
  Utmp(C<cmin) = -inf;
  U0(:,:,i) = Utmp;
  C = C - k;
  Utmp = util(C);
  Utmp(C<cmin) = -inf;
  U1(:,:,i) = Utmp;
end
clear Utmp


%% SOLUTION

% Intialize Value and Policy Functions
v  = -inf(n,T1+1,T2);
v0 = -inf(n,T1+1,T2);
v1 = -inf(n,T1+1,T2);
x0 = -inf(n,T1+1,T2);
x1 = -inf(n,T1+1,T2);


%% Solve Bellman Equation

tic

% Terminal Period
fprintf ('Termial Period')
for i=1:T1+1
  Cret = r*(R*X+P*(i-1));
  V = U0(:,:,i) + delta*util(Cret)/(1-delta);
  V(Cret<cmin) = -inf;
  [v0(:,i,T2),ix] = max(V,[],2);
  x0(:,i,T2) = x(ix);
  v(:,i,T2) = v0(:,i,T2);
end

% Other Periods
for t=T2-1:-1:1
  fprintf ('Period %5i\n',t)
  if t<=T1
    for i=1:t
      % conditional on not using contraception
      inxt = i+1;
      vnext = interp1(w(:),v(:,inxt,t+1),Wnxt)*q;
      vnext = vnext(:,ones(1,n))';
      [v0(:,i,t),ix] = max(U0(:,:,i)+delta*vnext,[],2);
      x0(:,i,t) = x(ix);
      % conditional on using contraception
      inxt = i;
      vnext = interp1(w,v(:,inxt,t+1),Wnxt)*q;
      vnext = vnext(:,ones(1,n))';
      [v1(:,i,t),ix] = max(U1(:,:,i)+delta*vnext,[],2);
      x1(:,i,t) = x(ix);
      % update value function
      v(:,i,t) = max(v1(:,i,t),v0(:,i,t));
    end
  else
    for i=1:T1
      % no need to use contraception
      inxt = i;
      vnext = interp1(w,v(:,inxt,t+1),Wnxt)*q;
      vnext = vnext(:,ones(1,n))';
      [v0(:,i,t),ix] = max(U0(:,:,i)+delta*vnext,[],2);
      x0(:,i,t) = x(ix);
      % update value function
      v(:,i,t) = v0(:,i,t);
    end
  end
end

toc


%% SIMULATION
 
% Initialize History Arrays
nrep = 100;
isimhist = zeros(nrep,T2);
wsimhist = zeros(nrep,T2);

% Intialize
isim = ones(nrep,1);           % No children
wsim = y(discrand(nrep,q));    % No net savings, just income

% Simulate
for t=1:T2
  % Record current state
  isimhist(:,t) = isim;
  wsimhist(:,t) = wsim;
  % Generate random income
  if t<T2
    ynxt = y(discrand(nrep,q));
  else
    ynxt = zeros(nrep,1);
  end
  if t<=T1
    vv0 = -inf(nrep,1);
    vv1 = -inf(nrep,1);
    for i=1:t+1
      j = find(isim==i);
      vv0(j) = interp1(w,v0(:,i,t),wsim(j));
      vv1(j) = interp1(w,v1(:,i,t),wsim(j));
    end
    for i=1:t+1
      j = find(isim==i&vv1>vv0);
      % use contraception
      wsim(j) = R*interp1(w,x1(:,i,t),wsim(j))+ynxt(j);
      j = find(isim==i&vv0>=vv1);
      % do not use contraception
      wsim(j) = R*interp1(w,x0(:,i,t),wsim(j))+ynxt(j);
      isim(j) = isim(j)+1;
    end
  else
    % no need to use contraception
    for i=1:T1+1
      j = find(isim==i);
      wsim(j) = R*interp1(w,x0(:,i,t),wsim(j))+ynxt(j);
    end
  end
end
isimhist = [isimhist isim];
wsimhist = [wsimhist wsim];

% Plot Simulated Number of Children
figure
hold on
bar(isimhist(1:3,1:T1+1)'-1)
title('Number of Children at Start of Period')
xlabel('Period')
ylabel('Children')
xlim([1 T1+1])
set(gca,'ytick',[0:T1])

% Plot Simulated Wealths
figure
hold on
plot(1:T2+1,wsimhist(1:3,:),'LineWidth',3)
title('Wealth at Start of Period')
xlabel('Period')
ylabel('Wealth')
xlim([1 T2+1])

% Plot Simulated Number of Children
figure
hold on
bar(mean(isimhist(:,1:T1+1))-1)
title('Mean Number of Children at Start of Period')
xlabel('Period')
ylabel('Children')
xlim([1 T1+1])
set(gca,'ytick',[0:T1])

% Plot Simulated Wealths
figure
hold on
plot(1:T2+1,mean(wsimhist),'LineWidth',3)
title('Mean Wealth at Start of Period')
xlabel('Period')
ylabel('Wealth')
xlim([1 T2+1])
