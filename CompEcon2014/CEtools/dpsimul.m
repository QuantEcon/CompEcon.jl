%% DPSIMUL
%
%  Discrete-Time Single-Agent Dynamic Optimization Model Simulator
%
%  Performs Monte-Carlo simulation of the controlled Markov process 
%  of a discrete-time dynamic optimization model solved by dpsolve.
%
%  Usage
%    [ssim,xsim,isim,jsim] = dpsimul(model,basis,nper,sinit,iinit,s,v,x)
%  Let
%    ds   = dimension of continuous state s
%    dx   = dimension of continuous action x
%    de   = dimension of continuous state transition shock e
%    ni   = number of discrete states i
%    nj   = number of discrete actions j
%    ns   = number of continuous state nodes
%    nrep = number of replications simulated
%    nper = number of periods simulated
%  Input
%    model : structured array containing model specifications
%    basis : ds-dimensional basis for functions defined on continuous
%            state space
%    nper  : number of periods simulated, a positive integer
%    sinit : nrep.ds initial continuous states
%    iinit : nrep.1 initial discrete states
%    s     : ns.ds  continuous state nodes
%    v     : ns.ni.nj optimal values at continuous state nodes, per 
%            discrete state and action
%    x     : ns.dx.ni.nj optimal continuous actions at continuous state 
%            nodes, per discrete state and action
%  Output
%    ssim  : nrep.nper+1.ds simulated continuous states
%    xsim  : nrep.nper.dx simulated continuous actions
%    isim  : nrep.nper+1 simulated discrete states
%    jsim  : nrep.nper simulated discrete actions
%    if nrep=1, arrays are squeezed along this singleton dimension

%  Copyright(c) 1997-2014
%    Mario J. Miranda - miranda.4@osu.edu
%    Paul L. Fackler  - paul_fackler@ncsu.edu

function [ssim,xsim,isim,jsim] = dpsimul(model,basis,nper,sinit,iinit,s,v,x)

% Set model fields to default if nonexistent (see above)
if ~isfield(model,'horizon'),  model.horizon = inf;              end  
if ~isfield(model,'func'),     error('Missing Function File');   end  
if ~isfield(model,'params'),   error('Missing Parameter List');  end  
if ~isfield(model,'discount'), error('Missing Discount Factor'); end  
if ~isfield(model,'ds'),       model.ds = 1;                     end  
if ~isfield(model,'dx'),       model.dx = 1;                     end  
if ~isfield(model,'ni'),       model.ni = 1;                     end  
if ~isfield(model,'nj'),       model.nj = 1;                     end  
if ~isfield(model,'e'),        model.e  = 0;                     end  
if ~isfield(model,'w'),        model.w  = 1;                     end  
if ~isfield(model,'q'),        model.q  = [];                    end
if ~isfield(model,'h'),        model.h  = [];                    end
if ~isfield(model,'X'),        model.X  = zeros(1,model.dx);     end
model.ni = max(1,model.ni);
model.nj = max(1,model.nj);

% Unpack model structure
T      = model.horizon;  
func   = model.func;  
params = model.params;
ds     = model.ds;  
dx     = model.dx;  
ni     = model.ni;  
nj     = model.nj;  
e      = model.e;   
w      = model.w;  
q      = model.q;
h      = model.h;

% Equivalent transition probability matrix for deterministic discrete state
if ni==1; h = ones(nj,1); end
if isempty(q)
  q = zeros(ni,ni,nj);
  if isempty(h)
    if ni==nj
      disp('Neither q nor h specified; will default to h(j,i)=j.')
      disp('Please specify q or h if this is not correct.')
      disp(' ')
      h = (1:nj)'*ones(1,ni);
    else
      error('Either q or h must be specified.')
    end
  end
  for i=1:ni
    for j=1:nj
      q(i,h(j,i),j) = 1;
    end
  end
end

% Determine number of collocation nodes, variable dimensions, number of
% replications and periods simulated
n  = basis.n;               % number of continuous state variable collocation
                            % ... node coordinates by dimension
nc = prod(n);               % number of continuous state variable collocation nodes
ns = size(s,1);             % number of continuous state variables input
nrep = size(sinit,1);       % number of replications
nper = min(nper,T);         % number of periods simulated cannot exceed horizon

% Compute collocation matrix and initialize basis coefficient arrays
Phi = funbase(basis,s);
cx  = zeros(nc,dx,ni,nj);
cv  = zeros(nc,ni,nj);

% Reshape x and v to conform to format recognized by this routine
if nargin<8, x = []; end
if T==inf
  x = reshape(x,ns,dx,ni,nj);
  v = reshape(v,ns,ni,nj);
else
  x = reshape(x,ns,dx,ni,nj,T);
  v = reshape(v,ns,ni,nj,T+1);
end

% Compute basis coefficient arrays for infinite horizon model
if T==inf
  for i=1:ni
    for j=1:nj
      cx(:,:,i,j) = Phi\x(:,:,i,j);
      cv(:,i,j)   = Phi\v(:,i,j);
    end
  end
end

% Initialize output arrays
ssim = zeros(nrep,nper+1,ds);
xsim = zeros(nrep,nper,dx);
isim = zeros(nrep,nper+1);
jsim = zeros(nrep,nper);

% Initialize states
ss = sinit;
ii = iinit;
ssim(:,1,:) = ss;
isim(:,1)   = ii;

% Simulate the model
for ip=1:nper
  if T<inf
    for i=1:ni
      for j=1:nj
        cx(:,:,i,j) = Phi\x(:,:,i,j,ip);
        cv(:,i,j)   = Phi\v(:,i,j,ip);
      end
    end
  end
  xx = zeros(nrep,dx);
  jj = zeros(nrep,1);
  vv = zeros(nrep,nj);
  for i=1:ni
    ir = find(ii==i);
    if isempty(ir), continue, end
    vv(ir,:) = funeval(squeeze(cv(:,i,:)),basis,ss(ir,:));
  end
  [vm,jmax] = max(vv,[],2);
  for i=1:ni
    for j=1:nj
      ir = find(ii==i&jmax==j);
      if isempty(ir), continue, end
      if dx>0
        cc = squeeze(cx(:,:,i,j));
        xx(ir,:) = funeval(cc,basis,ss(ir,:));
        [xl,xu]  = feval(func,'b',ss(ir,:),xx(ir,:),i,j,[],[],params{:});
        xx(ir,:) = min(max(xx(ir,:),xl),xu);
      end
      jj(ir) = j;
    end
  end
  if length(w)==1
     ee = e(ones(nrep,1),:);
  else
     ee = e(discrand(nrep,w),:);
  end
  iiold = ii;
  ii = ones(nrep,1);
  for i=1:ni
    for j=1:nj
      ir = find(iiold==i&jj==j);
      if isempty(ir), continue, end
      if ni>1
        ii(ir) =  discrand(length(ir),q(i,:,j));
      end
      ss(ir,:) = feval(func,'g',ss(ir,:),xx(ir,:),i,j,ii(ir),ee(ir,:),model.params{:});
    end
  end
  ssim(:,ip+1,:) = ss;
  isim(:,ip+1)   = ii;
  xsim(:,ip,:)   = xx;
  jsim(:,ip)     = jj;
end

if T==inf
  ssim(:,nper+1,:) = [];
  isim(:,nper+1)   = [];
end

ssim = squeeze(ssim);
xsim = squeeze(xsim);
isim = squeeze(isim);
jsim = squeeze(jsim);