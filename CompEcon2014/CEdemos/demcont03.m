function demcont03

%% DEMCONT03 Optimal Economic Growth Model (Stochastic)

% Preliminary tasks
demosetup(mfilename)


%% FORMULATION

% Model Parameters
alpha =  0.4;       % capital share
delta =  0.05;      % capital depreciation rate
gamma =  1;         % technology mean reversion coefficient
theta =  2.0;       % relative risk aversion
rho   =  0.05;      % discount rate
sigma =  0.1 ;      % volatility


%% DETERMINISTIC STEADY-STATE

% Compute Deterministic Steady State
kst = ((delta+rho)/alpha)^(1/(alpha-1));
qst = kst^alpha-delta*kst;
pst = qst^-theta;
vst = ((qst.^(1-theta))/(1-theta))/rho;
xst = [kst;qst];
disp('Steady State')
disp(xst)


%% SOLVE BELLMAN EQUATION DIRECTLY USING BROYDEN

% Approximation Structure
n = [15 5];
smin = [ 6 0.9];
smax = [14 1.1];
basis = fundefn('cheb',n,smin,smax);
s = funnode(basis);

% Define Residual Function
k     = @(s) s(:,1);
y     = @(s) s(:,2);
V     = @(c,s) funeval(c,basis,s);
Vk    = @(c,s) funeval(c,basis,s,[1 0]);
Vy    = @(c,s) funeval(c,basis,s,[0 1]);
Vyy   = @(c,s) funeval(c,basis,s,[2 2]);
q     = @(c,s) Vk(c,s).^(-1/theta);
resid = @(c,s) (q(c,s).^(1-theta))/(1-theta) ...
       + (y(s).*k(s).^alpha-delta*k(s)-q(c,s)).*Vk(c,s) ...
       + gamma*(1-y(s)).*Vy(c,s) ...
       - 0.5*(sigma^2)*y(s).*Vyy(c,s) - rho*V(c,s);

% Solve Bellman Equation Directly Using BROYDEN
tic
v = (((rho*k(s)).^(1-theta))/(1-theta));
c = funfitxy(basis,s,v);
optset('broyden','showiters',1)
c = broyden(resid,c,s);
toc
norm(resid(c,s))

%% SOLVE BELLMAN EQUATION USING SCSOLVE

% Model Structure
model.func   = @func;
model.rho    = rho;
model.params = {alpha,delta,gamma,theta,sigma};

% Solve Bellman Equation Using SCSOLVE
tic
v = (((rho*k(s)).^(1-theta))/(1-theta));
[c,scoord,v,q,res] = scsolve(model,basis,v);
toc


%% PLOT SOLUTION

% Reduce to One Dimension
s = gridmake(scoord);
j = find(s(:,2)==1);
k = s(j,1);
v = v(j);
q = q(j);
res = res(j);

% Plot Value Function
figure
hold on
plot(k,v)
plot(kst,vst,'k*','LineWidth',8)
title('Value Function')
xlabel('Capital Stock')
ylabel('Value')

% Plot Shadow Price Function
figure
hold on
p = funeval(c,basis,s(j,:),[1 0]);
plot(k,p)
plot(kst,pst,'k*','LineWidth',8)
title('Shadow Price Function')
xlabel('Capital Stock')
ylabel('Shadow Price')

% Plot Optimal Policy
figure
hold on
plot(k,q)
plot(kst,qst,'k*','LineWidth',8)
title('Optimal Policy')
xlabel('Capital Stock')
ylabel('Consumption')

% Plot Residual
figure
plot(k,res,k,0*k,'k:')
title('Bellman Equation Residual')
xlabel('Capital Stock')
ylabel('Residual')

% 
% %% Save Plots as EPS Files
% printfigures(mfilename,4)


%% Function File for scsolve
function out = func(flag,s,q,Vs,alpha,delta,gamma,theta,sigma)
k = s(:,1);
y = s(:,2);
n = length(k);
switch flag
  case 'x'
    Vk = Vs(:,1);
    out = Vk.^(-1./theta);
  case 'f'
    out = (q.^(1-theta))./(1-theta);
  case 'g'
    out = [(y.*k.^alpha-delta*k-q)  gamma*(1-y)];
    % out = alpha*log(k+1)-delta*k-q;
  case 'sigma'
    out = zeros(n,2,2);
    out(:,2,2) = sigma*sqrt(y);
end