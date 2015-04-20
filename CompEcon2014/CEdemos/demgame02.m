function demgame02

%% DEMGAME02 Income Redistribution Game
%
% Two infinitely-lived agents who agree to share income risk must make 
% consumption and investment decisions. 
%
%    States
%        s       wealth, post income transfer
%    Actions
%        x       savings
%    Parameters
%        alpha   utility parameter
%        beta    production elasticity
%        gamma   capital survival rate
%        psi     risk sharing rate
%        sigma   production shock volatility
%        delta   discount factor

% Preliminary tasks
demosetup(mfilename)


%% Formulation

% Model Parameters
alpha = [0.3 0.3];                      % utility parameter
beta  = [0.5 0.5];                      % production elasticity
gamma = [0.9 0.9];                      % capital survival rate
psi   = 0.1;                            % risk sharing rate
sigma = [0.2 0.2];                      % production shock volatility
delta = 0.9;                            % discount factor

% Shock Distribution
nshocks = [3 3];                                % number of shocks
[e,w] = qnwlogn(nshocks,[0 0],diag(sigma.^2));  % shocks and probabilities

% Model Structure
model.func = @func;                     % model functions
model.params = {alpha beta gamma psi};  % other parameters
model.discount = delta;                 % discount factor
model.dp = 2;                           % number of players
model.ds = 2;                           % dimension of state variable s
model.dx = 1;                           % dimension of individual action variable x_i
model.e = e;                            % shocks
model.w = w;                            % probabilities

% Compute deterministic steady=state
xstar = ((1/delta-gamma)./beta).^(1./(beta-1));
sstar = gamma.*xstar+xstar.^beta;
ratio = xstar./sstar;
vstar = (1/(1-delta))*((sstar-xstar).^(1-alpha))./(1-alpha);
pstar = (sstar-xstar).^(-alpha);

% Approximation Structure
n      = [9 9];                         % degree of approximation
smin   = [ 2  2];                       % minimum state
smax   = [12 12];                       % maximum state
basis  = fundefn('cheb',n,smin,smax);   % basis functions
snodes = funnode(basis);                % state collocaton nodes

% Initialize arrays
xinit = [snodes(:,1)*ratio(1) snodes(:,2)*ratio(2)];
vinit = [vstar(1)+pstar(1)*(snodes(:,1)-sstar(1)) vstar(2)+pstar(2)*(snodes(:,2)-sstar(2))] ;

% Check Model Derivatives
gamecheck(model,(smax+smin)/2,ratio.*(smax+smin)/2);


%% Solution

% Solve Bellman Equation
optset('gamesolve','nres',3);
[c,s,v,x,resid] = gamesolve(model,basis,vinit,xinit);


%% Analysis

% Reshape Output for Plotting
nr = n*3;
s1 = reshape(s(:,1),nr);
s2 = reshape(s(:,2),nr);
v = reshape(v,[nr 2]);
x = reshape(x,[nr 2]);
resid = reshape(resid,[nr 2]);

% Compute Shadow Prices
p = zeros(prod(nr),2);
p(:,1) = funeval(c(:,1),basis,s,[1 0]);
p(:,2) = funeval(c(:,2),basis,s,[0 1]);
p = reshape(p,[nr 2]);

% Plot Optimal Policy
for ip=1:2
  figure
  surf(s1,s2,x(:,:,ip),'FaceColor','interp','EdgeColor','interp');
  title(['Optimal Investment Policy: Player ' num2str(ip)])
  xlabel('Wealth 1')
  ylabel('Wealth 2')
  zlabel('Investment')
  xlim([0 15])
  ylim([0 15])
end

% Plot Value Function
for ip=1:2
  figure
  surf(s1,s2,v(:,:,ip),'FaceColor','interp','EdgeColor','interp');
  title(['Value Function: Player ' num2str(ip)])
  xlabel('Wealth 1')
  ylabel('Wealth 2')
  zlabel('Value')
  xlim([0 15])
  ylim([0 15])
end

% Plot Own Shadow Price Function
for ip=1:2
  figure
  surf(s1,s2,p(:,:,ip),'FaceColor','interp','EdgeColor','interp');
  title(['Own Shadow Price: Player ' num2str(ip)])
  xlabel('Wealth 1')
  ylabel('Wealth 2')
  zlabel('Price')
  xlim([0 15])
  ylim([0 15])
end

% Plot Residual
for ip=1:2
  figure
  surf(s1,s2,resid(:,:,ip),'FaceColor','interp','EdgeColor','interp');
  title(['Approximation Residual: Player ' num2str(ip)])
  xlabel('Wealth 1')
  ylabel('Wealth 2')
  zlabel('Residual')
  xlim([0 15])
  ylim([0 15])
end

% Simulate Model
ss{1} = s1(:,1);
ss{2} = s2(1,:)';
nyrs = 20;
nrep = 2000;
sinit = 7.49*ones(nrep,2);
[spath,xpath] = gamesimul(model,sinit,nyrs,ss,x);

% Plot Expected State Path
for ip=1:2
  figure
  plot(0:nyrs,mean(squeeze(spath(:,ip,:))))
  title(['Expected Wealth: Player ' num2str(ip)])
  xlabel('Year')
  ylabel('Wealth')
end

% Plot Expected Policy Path
for ip=1:2
  figure
  plot(0:nyrs,mean(squeeze(xpath(:,ip,:))))
  title(['Expected Investment: Player ' num2str(ip)])
  xlabel('Year')
  ylabel('Investment')
end

% Save Plots as EPS Files
printfigures(mfilename,12)


%% Function File
%
%  Called by gamesolve
%
%  A user-supplied function that returns the bound, reward, and state 
%  transitions and their first and second derivatives with respect to 
%  player i's action x_i, at an arbitrary number ns of states s and joint 
%  actions x, according to the format
%    [out1,out2,out3] = func(flag,i,s,x,e,params)
%  Let
%    ns = number of state nodes
%    dp = number of players
%    ds = dimension of state s
%    dx = dimension of individual action x_i
%    de = dimension of shock e
%  Function File Input
%    flag      : flag indicating function to be evaluated 
%    i         : player index, an integer between 1 and dp
%    s         : ns.ds    states
%    x         : ns.dx.dp actions
%    e         : ns.de    shocks
%    params    : parameters passed to function file
%  Function File Output
%    if flag = 'b', returns bounds on individual action x_i
%      out1    : ns.dx lower bounds on individual action x_i
%      out2    : ns.dx upper bounds on individual action x_i
%      out3    : empty
%    if flag = 'f', returns reward and derivatives w.r.t. x_i
%      out1    : ns.1     reward f_i
%      out2    : ns.dx    first derivative of f_i with respect to x_i
%      out3    : ns.dx.dx second derivative of f_i with respect to x_i
%    if flag = 'g', returns state transition and derivatives w.r.t. x_i
%      out1    : ns.ds       state g
%      out2    : ns.ds.dx    first derivative of g with respect to x_i
%      out3    : ns.ds.dx.dx second derivative of g with respect to x_i

function [out1,out2,out3] = func(flag,i,s,x,e,alpha,beta,gamma,psi)

ns = size(s,1);
if ns>1
  x = squeeze(x);
end
switch flag
  case 'b'
    xl = zeros(ns,1);
    xu = 0.99*s(:,i);
    out1=xl; out2=xu; out3=[];
  case 'f'
    f   = ((s(:,i)-x(:,i)).^(1-alpha(i)))/(1-alpha(i));
    fx  = -(s(:,i)-x(:,i)).^(-alpha(i));
    fxx = -alpha(i)*(s(:,i)-x(:,i)).^(-alpha(i)-1);
    out1=f; out2=fx; out3=fxx;
  case 'g'
    g1   = gamma(1)*x(:,1) + (1-psi)*e(:,1).*x(:,1).^beta(1) + psi*e(:,2).*x(:,2).^beta(2);
    g2   = gamma(2)*x(:,2) + (1-psi)*e(:,2).*x(:,2).^beta(2) + psi*e(:,1).*x(:,1).^beta(1);
    if i==1
      dg1  = gamma(1) + beta(1)*(1-psi)*e(:,1).*x(:,1).^(beta(1)-1);
      dg2  = beta(1)*psi*e(:,1).*x(:,1).^(beta(1)-1);
      ddg1 = (beta(1)-1)*beta(1)*(1-psi)*e(:,1).*x(:,1).^(beta(1)-2);
      ddg2 = (beta(1)-1)*beta(1)*psi*e(:,1).*x(:,1).^(beta(1)-2);
    else
      dg1  = beta(2)*psi*e(:,2).*x(:,2).^(beta(2)-1);
      dg2  = gamma(2) + beta(2)*(1-psi)*e(:,2).*x(:,2).^(beta(2)-1);
      ddg1 = (beta(2)-1)*beta(2)*psi*e(:,2).*x(:,2).^(beta(2)-2);
      ddg2 = (beta(2)-1)*beta(2)*(1-psi)*e(:,2).*x(:,2).^(beta(2)-2);
    end
    g    = [g1   g2];
    gx   = [dg1  dg2];
    gxx  = [ddg1 ddg2];
    gx  = reshape(gx,ns,2,1);
    gxx = reshape(gxx,ns,2,1,1);
    out1=g; out2=gx; out3=gxx;
end
%%