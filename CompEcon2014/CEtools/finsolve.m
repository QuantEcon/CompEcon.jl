% FINSOLVE Solves continuous time asset pricing problems
% USAGE
%   [c,V]=finsolve(model,fspace,alg,s,N);
% INPUTS
%   model   : a model structure (see below)
%   fspace  : name of projection space structure
%   alg     : algorithm used (lines, explicit, implicit, CN, or stiff)
%   s       : collocation nodal values of the state
%   N       : number of time steps
% OUTPUTS
%   c       : value function projection coefficients n x 1 
%                (or n x N+1 with keepall option)
%   V       : asset value n x 1
%                (or n x N+1 with keepall option)
%
% The model structure has the following fields:
%    model.func     - the name of the model function file
%    model.T        - the time to maturity of the asset
%    model.american :  1 if there is an early exercise feature for buyer
%                     -1 if there is an early exercise feature for seller
%    model.barrier  - 2x(2+d) matrix: rows specify lower and upper barriers
%    model.params   - additional parameters to be passed to model.func
% To price barrier options model.barrier should have the form
%   [R_a alpha_a beta_a;
%    R_b alpha_b beta_b];
% where R_a is paid anytime beta_a*S<=alpha_a and 
%       R_b is paid anytime beta_b*S<=alpha_b.
%
% The model function file must have the following syntax:
%    function out=func(flag,S,additional parameters);
%      switch flag
%       case 'rho'
%         out= instantaneous risk-free interest rate
%       case 'mu'
%         out= drift on the state process
%       case 'sigma'
%         out = volatility on the state process
%       case 'delta'
%         out = the payout flow (dividend) on the derivative asset
%       case 'V0'
%         out = exercise value of the asset
%       end
%
% Note: The text documentation has an argument, t, which is not
%  incorporated in the current version and is not used in the examples.
%
% USER OPTIONS (SET WITH OPSET)
%   keepall   : keeps coefficient values (c) for all time steps
%                 (this may be limited by available memory)
%   payatT    : for barrier options, 1 if rebate paid at expiration

% Revision 9/1/2010 replaced the use of lusolve with the pair luget/lusol

% Copyright (c) 1997-2010, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function [c,V]=finsolve(model,fspace,alg,s,N)

keepall  = optget('finsolve','keepall',0);
payatT   = optget('finsolve','payatT',0);

if ~exist('alg','var') || isempty(alg), alg = 'lines';         end
if ~exist('s','var')   || isempty(s),   s   = funnode(fspace); end
if ~exist('N','var')   || isempty(N),   N   = 1;               end

mffile = model.func;
T      = model.T;
params = model.params;

n = prod(fspace.n);
d = fspace.d;
S = gridmake(s);         % expand the grid

if isfield(model,'american'), american=model.american;
else                          american=0;
end

% Set up barier problem, if appropriate
if isfield(model,'barrier')
  barrier=1;
  temp=model.barrier;
  Ra=temp(1,1); Rb=temp(2,1);        % barrier rewards
  Binda=S*temp(1,3:end)'<=temp(1,2); % points less than lower barrier
  Bindb=S*temp(2,3:end)'>=temp(2,2); % points greater than upper barrier
else
  barrier=0;
end

% Compute collocation matrix
bases=funbasex(fspace,s,[0;1;2]);
B=ctbasemake(feval(mffile,'mu',S,params{:}),bases,1);
B=B+ctbasemake(feval(mffile,'sigma',S,params{:}),bases,2);
rho=ctbasemake(feval(mffile,'rho',S,params{:}),bases,0);
B=B-rho;
if barrier
  if payatT 
    B(Binda,:)=-rho(Binda,:);
    B(Bindb,:)=-rho(Bindb,:);
  else
    B(Binda,:)=0;
    B(Bindb,:)=0;
  end
end
clear rho 
Phi=bases.vals(1,d:-1:1);
bases.vals([2;3],:)=[];

% Get dividend information
b=feval(mffile,'delta',S,params{:});
hasdiv=~isempty(b);

% Get initial (terminal) value and compute initial coefficients
V0=feval(mffile,'V0',S,params{:});
if barrier
  V0(Binda)=Ra;
  V0(Bindb)=Rb;
end
c0=ckronxi(Phi,V0);

% Define matrix operators (or call ODE solver for stiff method)
Delta=T/N;
mass=1;
tensor=0;
switch alg
case 'lines'
  B=ckronxi(Phi,B);                   % inv(Phi)*B
  A=expm(full(Delta*B));
  if hasdiv, a=(A-speye(n))*(B\ckronxi(Phi,b)); end
  mass=0;
case 'explicit'
  A=ckron(Phi)+Delta*B;
  if hasdiv, a=Delta*b;end
  tensor=1;
  mass=0;
case 'implicit'
  A=ckron(Phi);
  M=A-Delta*B;
  if hasdiv, a=Delta*b;end
  MLU=luget(M);
case {'CN','trapezoid'}
  B=(Delta/2)*B;
  A=ckron(Phi);
  M=A-B;
  A=A+B;
  if hasdiv, a=Delta*b; end
  MLU=luget(M);
case {'stiff','ode'}
  if alg(1)=='o'
    solver=@ode45;
  else
    solver=@ode15s;
  end
  Phi=ckron(Phi);
  if ~hasdiv, b=[]; end
  options=odeset('Jacobian','on','Jconstant','on','Mass','M');
  if keepall
    t=linspace(0,T,N+1)';
    [t,c]=solver('finode',t,c0,options,b,B,Phi);
    c=c';
  else
    [t,c]=solver('finode',[0 T],c0,options,b,B,Phi);
    c=c(end,:)';
  end
  V=Phi*c;
  return
otherwise
  error('Method option is invalid')
end
clear B M

% Invert the bases matrices to avoid repeated inversions
if tensor || american || barrier
  for i=1:d
    Phi{i}=inv(Phi{i});
    % if the inverse matrix is not sparse, convert it to full
    if nnz(Phi{i})/numel(Phi{i})>0.3,
      Phi{i}=full(Phi{i});
    end
  end
end

% Initialize coefficient arrays
c=c0;
if keepall
  C=zeros(n,N+1);
  C(:,1)=c0;
end

% Iterate over time
for i=2:N+1
  if hasdiv, c=a+A*c;
  else       c=A*c;
  end
  if mass
    if tensor, c=ckronx(Phi,c); 
    else       c=lusol(MLU,c); 
    end
  end
  if barrier || american
    V=funeval(c,fspace,bases);
    if barrier
      V(Binda)=Ra;
      V(Bindb)=Rb;
    end
    if american
      if american<0
        V=min(V0,V);
      else
        V=max(V0,V);
      end
    end
    c=ckronx(Phi,V);
  end
  if keepall, C(:,i)=c; end
end

if keepall, c=C; end
V=funeval(c,fspace,bases);
