
% CTSSRES Residual function for continuous time steady state solver
% See also ctsteadystate

% Copyright (c) 1997-2002, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function r=ctssres(s,model,basis,cv)
x=feval(model.func,'x',s,[],funeval(cv,basis,s,1),model.params{:});
r=feval(model.func,'g',s,x,[],model.params{:});