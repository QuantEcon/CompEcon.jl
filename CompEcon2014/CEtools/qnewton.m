%% QNEWTON 
%
%  Unconstrained finite-dimensional maximization
%
%  QNEWTON performs unconstrained finite-dimensional maximization using
%  choice of search direction and line search method.
%
%  Usage
%   [x,A] = qnewton(func,x,A,varargin)
%  Let
%    n          : dimension of domain
%  Input
%    func       : objective function of form fval=f(x,varargin)
%    x          : n.1 initial guess for local maximum
%    A          : n.n initial guess for inverse Hessian (optional)
%    varargin   : additional arguments for func (optional)
%  Output
%    x          : n.1 local maximum of func, if found
%    A          : n.n inverse Hessian approximation
%  Options
%    Options may be set by user with OPTSET (defaults in parentheses):
%    SearchMeth : 1)Steepest Ascent
%                 2)DFP 
%                 3)BFGS [default]
%    StepMeth   : 1)No search 
%                 2)STEPBHHH 
%                 3)STEPBT [default]
%                 4)Golden Search
%    maxit      : maximum major iterations [250]
%    maxstep    : maximum step search iterations [50]
%    tol        : convergence tolerence [sqrt(eps)]
%    eps0       : zero factor (used in convergence criteria) [1]
%    eps1       : zero factor (used in climbing criteria) [1e-12]
%    ShowIters  : 1 to show results at each iteration [0]

% Copyright (c) 1997-2014, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function [x,A] = qnewton(func,x,A,varargin)

SearchMeth = optget('qnewton','SearchMeth',3);
StepMeth   = optget('qnewton','StepMeth',3);
maxit      = optget('qnewton','maxit',250);
maxstep    = optget('qnewton','maxstep',50);
tol        = optget('qnewton','tol',sqrt(eps));
eps0       = optget('qnewton','eps0',1);
eps1       = optget('qnewton','eps1',1e-12);
ShowIters  = optget('qnewton','ShowIters',0);

n = size(x,1);
reset = 0;

f0 = feval(func,x,varargin{:});
g0 = fdjac(func,x,varargin{:});
if abs(g0)<eps, return; end

if nargin<3 || isempty(A)
    A = -eye(n)./max(abs(f0),1);
    reset = 1;
end

for it=1:maxit
    d = -(A*g0');                        % search direction
    if ((d'*g0')/(d'*d))<eps1            % must go uphill
        A = -eye(n)./max(abs(f0),1);     % otherwise use
        d = g0'./max(abs(f0),1);         % steepest ascent
        reset = 1;
    end
    [s,f] = optstep(StepMeth,func,x,f0,g0,d,maxstep,varargin{:});
    if f<=f0                             % Step search failure
        if reset
            warning('Iterations stuck in qnewton'), return
        else                             % Use steepest ascent
            A = -eye(n)./max(abs(f0),1);
            d = g0'./max(abs(f0),1);
            [s,f,err] =  optstep(StepMeth,func,x,f0,g0,d,maxstep,varargin{:});
            if err, warning('Cannot find suitable step in qnewton'), return; end
        end
    end
    d = s*d;
    x = x+d;
    if any(isnan(x)|isinf(x))
        error('NaNs or INFs encountered')
    end
    f = feval(func,x,varargin{:});
    g = fdjac(func,x,varargin{:});
    if ShowIters
        fprintf('qnewton: %4i %16.4f %16.4f %16.4f\n',it,f,norm(d),norm(g))
    end
    % Test convergence using Marquardt's criteria and gradient test
    if ((f-f0)/(abs(f)+eps0)<tol && all(abs(d)./(abs(x)+eps0)<tol)) ...
            || all(abs(g)<eps); return; end
    % Update inverse Hessian
    u = g'-g0'; ud = u'*d;
    if SearchMeth==1 || abs(ud)<eps      	% Steepest ascent
        A = -eye(n)./max(abs(f),1);
        reset = 1;
    elseif SearchMeth==2;                   % DFP update
        v = A*u;
        A = A + d*d'./ud - v*v'./(u'*v);
        reset = 0;
    elseif SearchMeth==3;                   % BFGS update
        w = d-A*u; wd = w*d';
        A = A + ((wd+wd') - ((u'*w)*(d*d'))./ud)./ud;
        reset = 0;
    end
    %  Update iteration
    f0 = f; 
    g0 = g;
end
warning('Maximum iterations exceeded in qnewton')


%% OPTSTEP
%
%  Compute step length for unconstrained finite-dimensional optimization
%
%  Usage
%   [s,f,errcode,iter] = optstep(method,func,x0,f0,g0,d,MaxIters,varargin)
%  Let
%    n          : dimension of domain
%  Input
%    func       : objective function being maximized
%    x0         : starting value
%    f0         : value of func at x0
%    g0         : gradient of func at x0
%    d          : search direction
%    MaxIter    : max itereations before trying something else
%    varargin   : additional arguments for func
%  Output
%    s          : optimal step in direction d
%    f          : value of func at x+s*d,
%    iter       : number of iterations used
%    errcode    : 0 if suitable step found

function [s,f,errcode,iter] = optstep(method,func,x0,f0,g0,d,MaxIters,varargin)

if nargin<6, method = 3; end
if nargin<7, MaxIters = 100; end

g0 = g0';
switch method
    case 1
        f = feval(func,x0+d,varargin{:});
        if f<f0
            s = 1; iter = 1; errcode = 0;
        else
            [s,f,iter,errcode] = stepgold(func,x0,f0,g0,d,MaxIters,varargin{:});
        end
    case 2
        [s,f,iter,errcode] = stepbhhh(func,x0,f0,g0,d,MaxIters,varargin{:});
        if errcode
            [s,f,iter2,errcode] = stepgold(func,x0,f0,g0,d,MaxIters,varargin{:});
            iter = iter+iter2;
        end
    case 3
        [s,f,iter,errcode] = stepbt(func,x0,f0,g0,d,MaxIters,varargin{:});
        if errcode
            [s,f,iter2,errcode] = stepgold(func,x0,f0,g0,d,MaxIters,varargin{:});
            iter = iter+iter2;
        end
    otherwise
        [s,f,iter,errcode] = stepgold(func,x0,f0,g0,d,MaxIters,varargin{:});
end


%% STEPBHHH
%
%  Computes an approximate minimum step length using algorithm discussed in
%  Berndt, et. al., Annals of Economic and Social Measurement, 1974, pp.
%  653-665.
%
%  Usage
%   [s,fs,iter,errcode] = stepbhhh(func,x0,f0,g0,d,MaxIters,varargin)
%  Input
%    func       : objective function being maximized
%    x0         : current value of x
%    f0         : value of func at x0
%    g0         : gradient of func at x0
%    d          : search direction       
%    MaxIters   : maximum function evaluations allowed
%    varargin   : additional arguments for func
%  Output
%    s          : optimal step in direction d
%    fs         : value of func at x+s*d,
%    iter       : number of iterations used
%    errcode    : 1 if maximum iterations exceeded

function [s,fs,iter,errcode] = stepbhhh(func,x0,f0,g0,d,MaxIters,varargin)

% Intializations
if nargin<6 || isempty(MaxIters), MaxIters = 25; end
delta = optget('optstep','bhhhcone',0.0001);
dg   = -d'*g0;                   % directional derivative
tol1 = dg*delta;
tol0 = dg*(1-delta);
s  = 1;
ds = 1;
errcode = 0;

% Bracket the cone
for iter=1:MaxIters
    x  = x0+s*d; 
    fs = feval(func,x,varargin{:});
    temp = (f0-fs)/s;
    if temp<tol0
        ds = 2*ds;
        s  = s+ds;
    else break
    end
end
if tol0<=temp && temp<=tol1, return, end

ds = ds/2;
s  = s-ds;
it = iter+1;

% Then use bisection to get inside it
for iter=it:MaxIters
    ds = ds/2;
    x  = x0+s*d;
    fs = feval(func,x,varargin{:});
    temp = (f0-fs)/s;
    if     temp>tol1, s = s-ds;
    elseif temp<tol0, s = s+ds;
    else return
    end
end
errcode = 1;


%% STEPBT
%
%  Computes an approximate minimum step length using algorithm similar to
%  Algorithm 6.3.5 in Dennis and Schnabel, Numerical Methods for 
%  Unconstrained Optimization and Nonlinear Equations or LNSRCH in sec 9.7
%  of Press, et al., Numerical Recipes.
%
%  Usage
%    [s,fs,iter,errcode] = stepbt(func,x0,f0,g0,d,MaxIters,varargin)
%  Input
%    func       : objective function being maximized
%    x0         : current value of x
%    f0         : value of func at x0
%    g0         : gradient of func at x0
%    d          : search direction       
%    MaxIters   : maximum backsteps allowed
%    varargin   : additional arguments for func
%  Output
%    s          : optimal step in direction d
%    fs         : value of func at x+s*d,
%    iter       : number of iterations used
%    errcode    : 1 if maximum iterations exceeded
%                 2 if cubic approximation finds negative root

function [s,fs,iter,errcode] = stepbt(func,x0,f0,g0,d,MaxIters,varargin)

% Intializations
delta = 1e-4;    % Defines cone of convergence; must be on (0,1/2)
ub = 0.5;        % Upper bound on acceptable reduction in s.
lb = 0.1;        % Lower bound on acceptable reduction in s.

errcode = 0;
dg = -d'*g0;                      % directional derivative

tol1 = delta*dg;
tol0 = (1-delta)*dg;

% full step
s = 1;
fs = feval(func,x0+d,varargin{:});
if -fs+f0 <= tol1, iter = 1; return, end

% quadratic approximation
s2 = s; fs2 = fs;
s  = -0.5*dg./(-fs+f0-dg);
s  = max(s,lb);
fs = feval(func,x0+s*d,varargin{:});
temp = (-fs+f0)/s;
if tol0<=temp && temp <= tol1, iter = 2; return, end

% cubic approximation
for iter=3:MaxIters
    temp = (s-s2)*[s*s;s2*s2];
    temp = [(-fs+f0-dg*s);(-fs2+f0-dg*s2)]./temp;
    a    = temp(1)-temp(2);
    b    = s*temp(2)-s2*temp(1);
    s2   = s; 
    fs2  = fs;
    if a==0                                 % quadratic fits exactly
        s = -0.5*dg/b;
    else
        disc = b*b - 3*a*dg;
        if disc<0, errcode = 2; return, end % complex root
        s = (sqrt(disc)-b)/(3*a);
    end
    s    = max(min(s,ub*s2),lb*s2);         % ensures acceptable step size
    fs   = feval(func,x0+s*d,varargin{:});
    temp = (-fs+f0)/s;
    if tol0<=temp && temp <= tol1, return, end
end

errcode = 1;


%% STEPGOLD
%
%  Computes an approximate minimum step length using golden search method.
%
%  Usage
%    [s,fs,iter,errcode] = stepgold(func,x0,f0,g0,d,MaxIters,varargin)
%  Input
%    func       : objective function being maximized
%    x0         : current value of x
%    f0         : value of func at x0
%    g0         : gradient of func at x0 (unused)
%    d          : search direction       
%    MaxIters   : maximum backsteps allowed
%    varargin   : additional arguments for func
%  Output
%    s          : optimal step in direction d
%    fs         : value of func at x+s*d,
%    iter       : number of iterations used
%    errcode    : 1 if maximum iterations exceeded

function [s,fs,iter,errcode] = stepgold(func,x0,f0,g0,d,MaxIters,varargin)

alpha1 = (3-sqrt(5))/2; alpha2 = (sqrt(5)-1)/2;

tol = 1e-4;                % tolerance used for Golden search algorithm
tol = tol*(alpha1*alpha2); % the bracket will be len/(alpha1*alpha2)
s = 1;
errcode = 1;               % 1 if the search is unsuccessful; otherwise 0
iter = 0;
s0 = 0;
%  Find a bracketing interval
fs = feval(func,x0+d,varargin{:});
if f0 >= fs, len = alpha1;
else
    for iter=1:MaxIters
        s = 2*s;
        fl = fs;
        fs = feval(func,x0+s*d,varargin{:});
        if fs <= fl; len = alpha1*(s-s0); break
        else s0 = s/2;
        end
    end
    if iter>=MaxIters, s = s/2; fs = fl; return; end
end

xl = x0+(s0+len)*d;
xs = x0+(s-len)*d;

s = s-len;
len = len*alpha2;  % len now measures relative distance between xl and xs

fs = feval(func,xs,varargin{:});
fl = feval(func,xl,varargin{:});
% Golden search to find minimum
while iter<MaxIters
    iter = iter+1;
    if fs<fl
        s = s-len;
        len = len*alpha2;
        xs = xl; xl = xl-len*d;
        fs = fl; fl = feval(func,xl,varargin{:});
    else
        len = len*alpha2;
        s = s+len;
        xl = xs; xs = xs+len*d;
        fl = fs; fs = feval(func,xs,varargin{:});
    end
    if len<tol, errcode = 0; break, end
end

if fl>fs, fs = fl; s = s-len; end