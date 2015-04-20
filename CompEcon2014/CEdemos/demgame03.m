function demgame03

%% DEMGAME03 Marketing Board Game
%
% The marketing boards of two countries that compete in a duopolistic 
% international agricultural commodity market must decide how much to 
% export and how much to store.
%
%    States
%        s       quantity available per country at beginning of period
%    Actions
%        x       stocks hend in inventory per country at end of period
%    Parameters
%        kappa   unit storage cost
%        gamma   inverse demand elasticity
%        xmax    maximum storage
%        mu      relative size of agents
%        yvol    yield volatility
%        delta   discount factor

% Preliminary tasks
demosetup(mfilename)


%% Formulation

% Model Parameters
kappa = 0.05;                           % unit storage cost
gamma = -0.5;                           % inverse demand elasticity
xmax  = [0.2 0.2];                      % maximum storage
mu    = log([0.5 0.5]);                 % relative size of agents
yvol  = 0.2;                            % yield volatility
delta = 0.95;                           % discount factor

% Shock Distribution
nshocks = [3 3];                        % number of shocks
cov     = (yvol^2)*eye(2);              % covariance matrix of log shocks
[e,w]   = qnwlogn(nshocks,mu,cov);      % shocks and proabilitities

% Model Structure
model.func = @func;                     % model functions
model.params = {kappa gamma xmax};      % other parameters
model.discount = delta;                 % discount factor
model.dp = 2;                           % number of players
model.ds = 2;                           % dimension of state variable s
model.dx = 1;                           % dimension of individual action variable x_i
model.e = e;                            % shocks
model.w = w;                            % probabilities

% Approximation Structure
n      = [20 20];                       % degree of approximation
smin   = min(e);                        % minimum supply
smax   = max(e)+xmax;                   % maximum supply
basis  = fundefn('spli',n,smin,smax);   % basis functions

% Check Model Derivatives
gamecheck(model,(smax+smin)/2,zeros(1,2));


%% Solution

% Solve Bellman Equation
optset('gamesolve','nres',3);
[c,s,v,x,resid] = gamesolve(model,basis);


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
  title(['Optimal Storage: Country ' num2str(ip)])
  xlabel('S_1');
  ylabel('S_2')
  zlabel('Storage')
  xlim([s1(1),s1(end)])
  ylim([s2(1),s2(end)])
  zlim([0 .12]);
end

% Plot Value Function
for ip=1:2
  figure
  surf(s1,s2,v(:,:,ip),'FaceColor','interp','EdgeColor','interp');
  title(['Value Function: Country ' num2str(ip)])
  xlabel('S_1')
  ylabel('S_2')
  zlabel('Value')
  xlim([s1(1),s1(end)])
  ylim([s2(1),s2(end)])
end

% Plot Own Shadow Price Function
for ip=1:2
  figure
  surf(s1,s2,p(:,:,ip),'FaceColor','interp','EdgeColor','interp');
  title(['Own Shadow Price: Country ' num2str(ip)])
  xlabel('S_1')
  ylabel('S_2')
  zlabel('Price')
  xlim([s1(1),s1(end)])
  ylim([s2(1),s2(end)])
end

% Plot Residual
for ip=1:2
  figure
  surf(s1,s2,resid(:,:,ip),'FaceColor','interp','EdgeColor','interp');
  title(['Own Shadow Price: Country ' num2str(ip)])
  xlabel('S_1')
  ylabel('S_2')
  zlabel('Residual')
  xlim([s1(1),s1(end)])
  ylim([s2(1),s2(end)])
end

% Simulate Model
ss{1} = s1(:,1);
ss{2} = s2(1,:)';
nyrs = 10;
nrep = 1000;
sinit = smax(ones(nrep,1),:);
[spath,xpath] = gamesimul(model,sinit,nyrs,ss,x);

% Plot Expected State Path
for ip=1:2
  figure
  plot(0:nyrs,mean(squeeze(spath(:,ip,:))))
  title(['Expected Supply: Country ' num2str(ip)])
  xlabel('Year')
  ylabel('Supply')
end

% Plot Expected Policy Path
for ip=1:2
  figure
  plot(0:nyrs,mean(squeeze(xpath(:,ip,:))))
  title(['Expected Storage: Country ' num2str(ip)])
  xlabel('Year')
  ylabel('Storage')
end

% Plot Market Price Function
figure
ss = gridmake(s);
xx = reshape(x,size(ss));
pp = sum(ss-xx,2).^gamma;
pp = reshape(pp,nr);
surf(s1,s2,pp,'FaceColor','interp','EdgeColor','interp');
title('Market Price')
xlabel('S_1');
ylabel('S_2')
zlabel('Market Price')
xlim([s1(1),s1(end)])
ylim([s2(1),s2(end)])
zlim([0.7 1.2])

% Save Plots as EPS Files
printfigures(mfilename,13)


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

function [out1,out2,out3] = func(flag,i,s,x,e,kappa,gamma,xmax)
ns = size(s,1);
if ns>1
  x = squeeze(x);
end
switch flag
  case 'b'
    xl = zeros(ns,1);
    xu = xmax(i)*ones(ns,1);
    out1=xl; out2=xu; out3=[];
  case 'f'
    q   = s-x;
    qtot = sum(q,2);
    p   = qtot.^gamma;
    px  = -gamma*qtot.^(gamma-1);
    pxx = (gamma-1)*gamma*qtot.^(gamma-2);
    f   = p.*q(:,i) - kappa*x(:,i);
    fx  = px.*q(:,i) - p - kappa;
    fxx = pxx.*q(:,i) - 2*px;
    out1=f; out2=fx; out3=fxx;
  case 'g'
    gx  = zeros(ns,2,1);
    gxx = zeros(ns,2,1,1);
    g   = x + e;
    gx(:,i) = ones(ns,1);
    out1=g; out2=gx; out3=gxx;
end
%%