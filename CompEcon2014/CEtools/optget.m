%% OPTGET
%
%  Utility to Retreive Current CompEcon Routine Options
%
%  Usage
%    optvalue = optget(funcname,optname,optvalue)
%  Input
%    funcname : name of function whose options are to be retrieved
%    optname  : name of option
%    optval   : value of option
%  Output
%    optval   : current value of the option
%
%  Called by CompEcon routines, not by user. If named field not
%    defined, it will be set to optvalue, but will not be changed
%    otherwise.
%
%  Use OPTSET %  to change a previously set field.
%
%  optget(funcname) returns current values of options structure.

%  Copyright(c) 1997-2014
%    Mario J. Miranda - miranda.4@osu.edu
%    Paul L. Fackler  - paul_fackler@ncsu.edu

function optvalue = optget(funcname,optname,optvalue)

funcname = lower(funcname);
optvar = [funcname '_options'];
eval(['global ' optvar])       % declare a global variable

if nargin==1                   % return the whole option structure
  optvalue = (eval(optvar));
  return
end

optname  = lower(optname);
% if structure empty or named field does not exist, set to passed value,
% otherwise return value in field
if isempty(eval(optvar)) | ~isfield(eval(optvar),optname)
  eval([optvar '.' optname '=optvalue;']);
else
  optvalue = eval([optvar '.' optname]);
end