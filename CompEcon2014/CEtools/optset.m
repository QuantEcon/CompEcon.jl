%% OPTSET
%
%  Utility for Overriding CompEcon Toolbox Routine Default Options
%
%  Usage
%    optset(funcname,optname,optvalue)
%  Input
%    funcname : name of function whose options are to be set by user
%    optname  : name of option
%    optval   : user-set value of option
%  Output
%    None. OPTSET operates on a global variable invisible to user.
%
% If optname='defaults', current option settings will be cleared. The next
%    time OPTSET is called, the default options will be restored.

%  Copyright(c) 1997-2014
%    Mario J. Miranda - miranda.4@osu.edu
%    Paul L. Fackler  - paul_fackler@ncsu.edu

function optset(funcname,optname,optvalue)

optvar = [lower(funcname) '_options'];      % name of global variable
optname  = lower(optname);                  % name of option field
if strcmp(optname,'defaults')
  eval(['clear global  ' optvar])           % clears global variable
else
  eval(['global  ' optvar])                 % declare global variable
  eval([optvar '.' optname '=optvalue;'])   % set specified field
end