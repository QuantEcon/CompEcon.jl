% FINODE ODE file for FINSOLVE

% Copyright (c) 1997-2001, Mario J. Miranda & Paul L. Fackler
% miranda.4@osu.edu, paul_fackler@ncsu.edu

function dc=finode(t,c,flag,b,B,Phi)
 switch flag
 case 'jacobian'
   dc=B;
 case 'mass'
   dc=Phi;
 otherwise
   if isempty(b), dc=B*c; else, dc=b+B*c; end
 end