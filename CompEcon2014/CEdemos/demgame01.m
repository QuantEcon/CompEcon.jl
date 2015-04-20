function demgame01

%% DEMGAME01 Production Capacity Game
%
% Two firms competing in a duopolistic product market must decide how much
% to invest in productive capacity.
%    States
%        s       capital stock at beginning of period
%    Actions
%        x       investent this period
%    Parameters
%        alpha   parameters
%        beta    parameters
%        gamma   parameters
%        psi     depreciation rate
%        delta   discount factor


% Preliminary tasks
demosetup(mfilename)


%% Formulation

% Model Parameters
alpha = [8.0 4.0];                      % parameters
beta  = [1.8 0.2];                      % parameters
gamma = [0.4 3.0];                      % parameters
psi   = 0.1;                            % depreciation rate
delta = 0.9;                            % discount factor

% Model Structure
model.func = @func;                     % model functions
model.discount = delta;                 % discount factor
model.params = {alpha,beta,gamma,psi};  % other parameters
model.dp = 2;                           % number of players
model.ds = 2;                           % dimension of state variable s
model.dx = 1;                           % dimension of individual action variable x_i

% Approximation Structure
n      = [8 8];                         % degree of approximation
smin   = [0.5 0.5];                     % minimum state
smax   = [1.5 1.5];                     % maximum state
basis  = fundefn('cheb',n,smin,smax);   % basis functions

% Check Model Derivatives
gamecheck(model,(smax+smin)/2,zeros(1,2));


%% Solution

% Solve Bellman Equation
optset('gamesolve','nres',4);
[c,s,v,x,resid] = gamesolve(model,basis);


%% Analysis

% Reshape Output for Plotting
nr = n*4;
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
  xlabel('S_1')
  ylabel('S_2')
  zlabel('Investment')
  zlim([0 .25])
end

% Plot Value Function
for ip=1:2
  figure
  surf(s1,s2,v(:,:,ip),'FaceColor','interp','EdgeColor','interp');
  title(['Value Function: Player ' num2str(ip)])
  xlabel('S_1')
  ylabel('S_2')
  zlabel('Value')
end

% Plot Own Shadow Price Function
for ip=1:2
  figure
  surf(s1,s2,p(:,:,ip),'FaceColor','interp','EdgeColor','interp');
  title(['Own Shadow Price: Player ' num2str(ip)])
  xlabel('S_1')
  ylabel('S_2')
  zlabel('Price')
end

% Plot Residual
for ip=1:2
  figure
  surf(s1,s2,resid(:,:,ip),'FaceColor','interp','EdgeColor','interp');
  title(['Approximation Residual: Player ' num2str(ip)])
  xlabel('S_1'); 
  ylabel('S_2')
  zlabel('Residual')
end

% Simulate Model
ss{1} = s1(:,1);
ss{2} = s2(1,:)';
nyrs = 30;
nrep =  2;
sinit = smin(ones(nrep,1),:);
[spath,xpath] = gamesimul(model,sinit,nyrs,ss,x);

% Plot Expected State Path
for ip=1:2
  figure
  plot(0:nyrs,mean(squeeze(spath(:,ip,:))))
  title(['Expected Capital Stock: Player ' num2str(ip)])
  xlabel('Year')
  ylabel('Capital')
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
    xu = inf*ones(ns,1);
    out1=xl; out2=xu; out3=[];
  case 'f'
    c = beta(1) + beta(2)./s;
    pi = ((alpha(1)-2*c(:,i)+c(:,3-i)).^2)/(9*alpha(2));
    f  = pi-(gamma(1)*x(:,i)+0.5*gamma(2)*x(:,i).^2);
    fx  = -(gamma(1)+gamma(2)*x(:,i));
    fxx = zeros(ns,1)-gamma(2);
    out1=f; out2=fx; out3=fxx;
  case 'g'
    gx  = zeros(ns,2,1);
    gxx = zeros(ns,2,1,1);
    g = (1-psi)*s + x;
    gx(:,i) = ones(ns,1);
    out1=g; out2=gx; out3=gxx;
end
%%