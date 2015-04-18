%% DPCHECK
%
%  Checks analytic derivatives supplied in function file passed to dpsolve.
%
%  Usage
%    dpcheck(model,s,x,i,j,in)
%  Let
%    ds   = dimension of continuous state variable s
%    dx   = dimension of continuous action variable x
%    de   = dimension of shock variable e
%    ni   = number of discrete states i
%    nj   = number of discrete actions j
%  Input
%    model : model structure
%    s     : 1.ds continuous state
%    x     : 1.dx continuous action
%    i     : discrete state, an integer between 1 and ni
%    j     : discrete action, an integer between 1 and nj
%    in    : discrete state next period, an integer between 1 and ni
%  See
%    dpsolve
 
%  Copyright(c) 1997-2011
%   Mario J. Miranda - miranda.4@osu.edu
%   Paul L. Fackler  - paul_fackler@ncsu.edu

function dpcheck(model,s,x,i,j,in)

if nargin<4, i =1; end;
if nargin<5, j =1; end;
if nargin<6, in=1; end;
if isfield(model,'w')
  e = model.w'*model.e; 
else
  e = 0;
end

[f,fx,fxx] = feval(model.func,'f',s,x,i,j,in,e,model.params{:});
[g,gx,gxx] = feval(model.func,'g',s,x,i,j,in,e,model.params{:});

fxfd  = fjac(model.func,3,'f',s,x,i,j,in,e,model.params{:});
gxfd  = fjac(model.func,3,'g',s,x,i,j,in,e,model.params{:});
fxxfd = fjac(model.func,[3 2],'f',s,x,i,j,in,e,model.params{:});
gxxfd = fjac(model.func,[3 2],'g',s,x,i,j,in,e,model.params{:});

err1 = max(abs(fx(:)-fxfd(:)));
err2 = max(abs(fxx(:)-fxxfd(:)));
err3 = max(abs(gx(:)-gxfd(:)));
err4 = max(abs(gxx(:)-gxxfd(:)));

if max([err1 err2 err3 err4])<1.e-4
   disp('Function File Passed Derivative Consistency Check')
else
   disp('Function File Failed Derivative Consistency Check')
end

if max(err1)>1.e-4
   disp('WARNING: Apparent Error in Derivative fx')
   disp('Max Discrepancy = ')
   disp(err1)
end

if max(err2)>1.e-4
   disp('WARNING: Apparent Error in Derivative fxx')
   disp('Max Discrepancy = ')
   disp(err2)
end

if max(err3)>1.e-4
   disp('WARNING: Apparent Error in Derivative gx')
   disp('Max Discrepancy = ')
   disp(err3)
end

if max(err4)>1.e-4
   disp('WARNING: Apparent Error in Derivative gxx')
   disp('Max Discrepancy = ')
   disp(err4)
end