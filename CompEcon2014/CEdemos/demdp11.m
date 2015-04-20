function demdp11

%% DEMDP11 Monetary Policy Model
%
% A central bank must set nominal interest rate so as to minimize
% deviations of inflation rate and GDP gap from established targets.
%
% States
%     s1      GDP gap
%     s2      inflation rate
% Actions
%     x       nominal interest rate
% Parameters
%     alpha   transition function constant coefficients
%     beta    transition function state coefficients
%     gamma   transition function action coefficients
%     omega   central banker's preference weights
%     sbar    equilibrium targets
%     delta   discount factor

% Preliminary tasks
demosetup(mfilename)


%% FORMULATION

% Model Parameters
alpha   = [0.9; -0.1];                  % transition function constant coefficients
beta    = [-0.5 0.2;0.3 -0.4];          % transition function state coefficients
gamma   = [-0.1; 0.0];                  % transition function action coefficients
omega   = [1 0;0 1];                    % central banker's preference weights
sbar    = [1;0];                        % equilibrium targets
sigma   = 0.08*eye(2);                  % shock covariance matrix
delta   = 0.9;                          % discount factor

% Compute Unconstrained Deterministic Steady-State
F0  = -0.5*sbar'*omega*sbar;
Fs  = sbar'*omega;
Fx  = 0;
Fss = -omega;
Fsx = [0 0]';
Fxx = 0;
G0  = alpha;
Gs  = beta;
Gx  = gamma;
[X,P,G,sstar,xstar] = lqsolve(F0,Fs,Fx,Fss,Fsx,Fxx,G0,Gs,Gx,delta);

% If Nonnegativity Constraint Violated, Re-Compute Deterministic Steady-State
if xstar<0
  I = eye(2,2);
  xstar = 0;
  sstar = (I-beta)\alpha;
%   lstar = (I-delta*beta')\(omega*(sstar-sbar));
%   gamma'*lstar
%   if gamma'*lstar>0, 'Something Wrong Here', end
end

% Reorient Deterministic Steady State
sstar   = sstar';
xstar   = xstar';
 
% Continuous State Shock Distribution
m   = [3 3];                            % number of shocks
mu  = [0 0];                            % means of shocks
[e,w] = qnwnorm(m,mu,sigma);            % shocks and probabilities

% Model Structure
model.func = @func;                     % model functions
model.params = {alpha beta gamma omega sbar};	% function parameters
model.discount = delta;                 % discount factor
model.ds = 2;                           % dimension of continuous state
model.dx = 1;                           % dimension of continuous action
model.ni = 0;                           % number of discrete states
model.nj = 0;                           % number of discrete actions
model.e  = e;                           % continuous state shocks
model.w  = w;                           % continuous state shock probabilities

% Approximation Structure
n = [30 30];                            % number of collocation coordinates, per dimension
smin = [-2 -3];                         % minimum states
smax = [ 2  3];                         % maximum states
basis = fundefn('spli',n,smin,smax);    % basis functions
s = funnode(basis);                     % collocation nodes

% Check Model Derivatives
dpcheck(model,sstar,xstar);


%% SOLUTION

% Solve Bellman Equation
optset('dpsolve','nr',5);
[c,s,v,x,resid] = dpsolve(model,basis);

% Reshape Output for Plotting
n = n*5;
s1 = reshape(s(:,1),n);
s2 = reshape(s(:,2),n);
v = reshape(v,n);
x = reshape(x,n);
resid = reshape(resid,n);

% Compute Shadow Prices
p1 = funeval(c,basis,s,[1 0]);
p2 = funeval(c,basis,s,[0 1]);
p1 = reshape(p1,n);
p2 = reshape(p2,n);

% Plot Optimal Policy
figure
surf(s1,s2,x,'FaceColor','interp','EdgeColor','interp')
title('Optimal Monetary Policy')
xlabel('GDP Gap')
ylabel('Inflation Rate')
zlabel('Nominal Interest Rate')
  
% Plot Value Function
figure
surf(s1,s2,v,'FaceColor','interp','EdgeColor','interp')
title('Value Function')
xlabel('GDP Gap');
ylabel('Inflation Rate')
zlabel('Value')
   
% Plot Shadow Price Function 1
figure
surf(s1,s2,p1,'FaceColor','interp','EdgeColor','interp')
title('Shadow Price of GDP Gap')
xlabel('GDP Gap')
ylabel('Inflation Rate')
zlabel('Price')
% Plot Shadow Price Function 2
figure
surf(s1,s2,p2,'FaceColor','interp','EdgeColor','interp')
title('Shadow Price of Inflation Rate')
xlabel('GDP Gap')
ylabel('Inflation Rate')
zlabel('Price')

% Plot Residual
figure
surf(s1,s2,resid,'FaceColor','interp','EdgeColor','interp')
title('Bellman Equation Residual')
xlabel('GDP Gap')
ylabel('Inflation Rate')
zlabel('Residual')


%% SIMULATION

% Simulate Model
rand('seed',0.945);
nper = 21;
nrep = 10000;
sinit = smax(ones(nrep,1),:);
iinit = ones(nrep,1);
[ssim,xsim] = dpsimul(model,basis,nper,sinit,iinit,s,v,x);
s1sim = ssim(:,:,1);
s2sim = ssim(:,:,2);

% Compute Ergodic Moments
s1avg = mean(s1sim(:)); 
s2avg = mean(s2sim(:)); 
xavg = mean(xsim(:)); 
s1std = std(s1sim(:)); 
s2std = std(s2sim(:)); 
xstd = std(xsim(:)); 

% Print Steady-State and Ergodic Moments
fprintf('Deterministic Steady-State\n') 
fprintf('   GDP Gap               = %5.2f\n'    ,sstar(1))
fprintf('   Inflation Rate        = %5.2f\n'    ,sstar(2))
fprintf('   Nominal Interest Rate = %5.2f\n\n'  ,xstar)
fprintf('Ergodic Means\n') 
fprintf('   GDP Gap               = %5.2f\n'    ,s1avg)
fprintf('   Inflation Rate        = %5.2f\n'    ,s2avg)
fprintf('   Nominal Interest Rate = %5.2f\n\n'  ,xavg)
fprintf('Ergodic Standard Deviations\n') 
fprintf('   GDP Gap               = %5.2f\n'    ,s1std)
fprintf('   Inflation Rate        = %5.2f\n'    ,s2std)
fprintf('   Nominal Interest Rate = %5.2f\n\n'  ,xstd)

% Plot Simulated and Expected State Paths 1
figure
hold on
plot(0:nper-1,s1sim(1:3,:))
plot(0:nper-1,mean(s1sim),'k')
plot(nper-1,s1avg,'k*')
title('Simulated and Expected State Paths')
xlabel('Period')
ylabel('GDP Gap')

% Plot Simulated and Expected State Paths 2
figure
hold on
plot(0:nper-1,s2sim(1:3,:))
plot(0:nper-1,mean(s2sim),'k')
plot(nper-1,s2avg,'k*')
title('Simulated and Expected State Paths')
xlabel('Period')
ylabel('Inflation Rate')

% Plot Simulated and Expected Policy Paths
figure
hold on
plot(0:nper-1,xsim(1:3,:))
plot(0:nper-1,mean(xsim),'k')
plot(nper-1,xavg,'k*')
title('Simulated and Expected Policy Paths')
xlabel('Period')
ylabel('Nominal Interest Rate')


%% Save Plots as EPS Files
printfigures(mfilename,8)


%% DPSOLVE FUNCTION FILE
%
%    User-supplied function called by dpsolve that returns the bound,
%    reward, and continuous state transition function values and
%    derivatives with respect to the continuous action x at an arbitrary
%    number ns of states and actions aaccording to the format
%      [out1,out2,out3] = func(flag,s,x,i,j,in,e,<params>)
%    For further information regarding input format, type > help dpsolve or
%    reference the DPSOLVE Handbook.

function [out1,out2,out3] = func(flag,s,x,i,j,in,e,alpha,beta,gamma,omega,sbar)

[n ds]  = size(s);
dx = 1;
switch flag
case 'b'
  out1  = -0*ones(n,1);
  out2  = inf*ones(n,1);
  out3 = [];
case 'f'
  s = s - sbar(:,ones(1,n))';
  out1 = zeros(n,1);
  for i=1:2
    for j=1:2
      out1 = out1 - 0.5*omega(i,j)*s(:,i).*s(:,j);
    end
  end
  out2 = zeros(n,1);
  out3 = zeros(n,1);  
case 'g'
  out1 = alpha(:,ones(1,n))' + s*beta' + x*gamma' + e;
  out2 = gamma(:,ones(1,n))';
  out3 = zeros(n,ds,dx,dx);
  out2 = reshape(out2,n,ds,dx);
end
%%