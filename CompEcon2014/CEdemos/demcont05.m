function demcont05

%% DEMCONT07 Production-Adjustment Model (Stochastic)
%
% Old demsc03

% Preliminary tasks
demosetup(mfilename)


%% FORMULATION

% Model Parameters
kappa = 0.5;          % speed of price mean reversion
alpha = 1;            % mean price
sigma = 0.2;          % price volatility
a     = 4;            % adjustment cost parameter
k     = 2;            % variable cost parameter
rho   = 0.1;          % discount rate

% Set Plotting Window
pmin = 0.4; 
pmax = 2.5;


%% SOLVE USING SCSOLVE

% Approximation Structure
n = [10 10];                           % number of state collocation coordinates
smin = [pmin 0];                       % minimum states
smax = [pmax 3];                       % maximum states
basis = fundefn('cheb',n,smin,smax);   % basis functions

% Model Structure
model.func   = @func;
model.rho    = rho;
model.params = {kappa alpha sigma a k};

% Solve Bellman Equation
[cv,s,v,x,resid] = scsolve(model,basis);
v = reshape(v,size(s{1},1),size(s{2},1));
x = reshape(x,size(s{1},1),size(s{2},1));
resid = reshape(resid,size(s{1},1),size(s{2},1));

% Compute Shadow Prices
n = [length(s{1}) length(s{2})];
p1 = funeval(cv,basis,s,[1 0]);
p2 = funeval(cv,basis,s,[0 1]);
p1 = reshape(p1,n);
p2 = reshape(p2,n);

% Plot Optimal Policy (Surface)
figure
hh = surf(s{1},s{2},x');
set(hh,'FaceColor','interp','EdgeColor','interp')
title('Optimal Production Policy')
xlabel('Price')
ylabel('Production Rate')
zlabel('Production Adjustment Rate')
view(-54,30)
xlim([.4 2.5])

% Plot Value Function (Surface)
figure
hh = surf(s{1},s{2},v');
set(hh,'FaceColor','interp','EdgeColor','interp')
title('Value Function')
xlabel('Price')
ylabel('Production Rate')
zlabel('Value')
view(-54,30)
xlim([.4 2.5])

% Plot Shadow Price Function 1 (Surface)
figure
hh = surf(s{1},s{2},p1');
set(hh,'FaceColor','interp','EdgeColor','interp')
title('Shadow Price of Market Price')
xlabel('Price')
ylabel('Production Rate')
zlabel('Shadow Price')
view(-54,30)
xlim([.4 2.5])

% Plot Shadow Price Function 2 (Surface)
figure
hh = surf(s{1},s{2},p2');
set(hh,'FaceColor','interp','EdgeColor','interp')
title('Shadow Price of Production Level')
xlabel('Price')
ylabel('Production Rate')
zlabel('Shadow Price')
view(-54,30)

% Plot Residual
figure
hh = surf(s{1},s{2},resid');
set(hh,'FaceColor','interp','EdgeColor','interp')
title('Approximation Residual')
xlabel('Price')
ylabel('Production Rate')
zlabel('Residual')
view(-54,30)
xlim([.4 2.5])


%% ERGODIC DISTRIBUTION

% Model Structure for Computing Ergodic Distribution of Price
model.func = @func2;                       % model functions
model.r = rho;                             % model functions
model.params = {kappa alpha sigma};        % other parameters

% Evaluate Ergodic Distribution
n = 50;                                    % number of state collocation coordinates
basis = fundefn('cheb',n,pmin,pmax);       % basis functions
cp = itodensity(model,basis);
p  = nodeunif(201,pmin,pmax);

% Plot Ergodic Distribution
figure
plot(p,funeval(cp,basis,p));
title('Ergodic Distribution of Price')
xlabel('Price')


%% Save Plots as EPS Files
printfigures(mfilename,6)


% Model function file for Production-Adjustment Model
function out1 = func(flag,s,x,Vs,kappa,alpha,sigma,a,k)
n = size(s,1);
switch flag
  case 'x';
    out1 = Vs(:,2)/a;
  case 'f';
    out1 = s(:,1).*s(:,2) - 0.5*k*s(:,2).^2 - 0.5*a*x.^2;
  case 'g';
    out1 = [kappa*(alpha-s(:,1)) x];
  case 'sigma'
    out1 = [sigma*s(:,1) zeros(n,3)];
end

% Supplemental model function file for Production-Adjustment Model
% Used to solve for stationary density of price
function out = func2(flag,s,kappa,alpha,sigma)
switch flag
  case 'mu';
    out = kappa*(alpha-s);
  case 'sigma'
    out = sigma*s;
end