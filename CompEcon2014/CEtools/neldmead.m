% NELDMEAD  Maximizes function via Nelder-Mead algorithm
% USAGE
%   [x,S] = neldmead(f,x,S,P1,P2,...);
% INPUTS
%   f         : string name of function fval=f(x)
%   x         : starting nx1 vector
%   S         : starting nx(n+1) simplex (optional)
%   P1,P2,... : additional parameters to pass to f
% OUTPUTS
%   x       : solution vector
%   S       : ending simplex sorted best to worst (x=S(:,1))
%
% Setable options (use OPTSET):
%   tol     : convergence tolerance
%   maxit   : maximum number of iterations

% Copyright (c) 1997-2002, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function [x,S] = neldmead(f,x,S,varargin)

maxit     = optget('neldmead','maxit',500);
tol1      = optget('neldmead','tol',sqrt(eps));
showiters = optget('neldmead','showiters',0);


n = length(x);
n1 = n+1;

% If a simplex is not passed, create one
if nargin<3 | isempty(S)
   S = [x x(:,ones(n,1))+diag(0.05*x+0.00025*(abs(x)<eps))];
end

% Evaluate at initial simplex
fS = zeros(1,n1);
for i=1:n1
  fS(i) = feval(f,S(:,i),varargin{:});
end
[temp,ind] = sort(fS);
rfactor = 1+2/n;

tol0=sqrt(eps);
tol2=tol1;

% Iterations
for iter=1:maxit
  besti = ind(end);
  worsti2 = ind(2);
  worsti = ind(1);
  if (fS(besti)-fS(worsti))/(abs(fS(worsti))+abs(fS(besti))+tol0)<tol1 & ...
    max(max(abs(S(:,besti+zeros(n1,1))-S),[],2)./(abs(S(:,besti))+tol0)) <= tol2 
    break
  end
  xworst=S(:,worsti);
  xtest = (2/n)*sum(S,2)-rfactor*xworst;   % get reflection
  fxtest = feval(f,xtest,varargin{:});
  if fxtest>fS(worsti2)                    % better than second worst
    if fxtest>fS(besti)                    % better than best
      xtest2 = 1.5*xtest-0.5*xworst;       % try expansion
      fxtest2 = feval(f,xtest2,varargin{:});
      if fxtest2>fxtest                    % accept expansion
        fxtest = fxtest2;
        xtest = xtest2;
      end
    end
    S(:,worsti) = xtest;
    fS(worsti) = fxtest;
  else                                     % reflection didn't change ordering
    xtest = 0.75*xworst+0.25*xtest;        % contract
    fxtest = feval(f,xtest,varargin{:});
    if fxtest>fS(worsti2)                  % accept contraction
      S(:,worsti) = xtest;
      fS(worsti) = fxtest;
    else                                   % shrink
      S = (S+S(:,besti+zeros(1,n1)))/2;
      for i=1:n1
        if i~=besti, fS(i)=feval(f,S(:,i),varargin{:}); end
      end
    end
  end
  [temp,ind] = sort(fS);
  if showiters, fprintf('%4i %6.2e %6.2e\n',[iter fS(1) fS(end)]); end
end

% Return in sorted order
ind=flipud(ind(:));
S=S(:,ind); x=S(:,1);

if iter>=maxit
  warning('Maximum iterations exceeded in Nelder-Mead')
end
