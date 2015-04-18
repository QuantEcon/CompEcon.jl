%% GAMECHECK 
%
%  Checks analytic derivatives supplied in function file passed to gamesolve.
% 
%  Usage
%    err = gamecheck(model,s,x,e)
%  Input
%    model : a game model structure
%    s     : a matrix of states
%    x     : a matrix of actions
%    e     : a matrix of shocks
%  Output
%    err   : a 4x1 vector containing maximal differences between analytic
%            and numerical derivatives df/dx, d^2f/dx^2, dg/dx, d^2g/dx^2, 
%            respectively.
%  See
%    gamesolve

%  Copyright(c) 1997-2010
%    Mario J. Miranda - miranda.4@osu.edu
%    Paul L. Fackler  - paul_fackler@ncsu.edu

function err=gamecheck(model,s,x,e)

if nargin<4;
  if isfield(model,'w')
    e=model.w'*model.e;
  else
    e=0;
  end
end

err = [];
for i=1:2
  [f,df,d2f] = feval(model.func,'f',i,s,x,e,model.params{:});
  [g,dg,d2g] = feval(model.func,'g',i,s,x,e,model.params{:});
  Jf = fjac(model.func,4,'f',i,s,x,e,model.params{:});
  Jg = fjac(model.func,4,'g',i,s,x,e,model.params{:});
  Hf = fjac(model.func,[4 2],'f',i,s,x,e,model.params{:});
  Hg = fjac(model.func,[4 2],'g',i,s,x,e,model.params{:});
  err1 = norm(df-Jf(i),inf);
  err2 = norm(d2f-Hf(i),inf);
  err3 = norm(dg-Jg(i,:),inf);
  err4 = norm(d2g-Hg(:,i)',inf);
  err  = [err err1 err2 err3 err4];
end

if max(err)>1.e-4
  disp('Possible Error in Derivatives')
  disp('Discrepancies in derivatives = ')
  disp(err)
end