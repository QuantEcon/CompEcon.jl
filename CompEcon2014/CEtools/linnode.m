% LINNODE Standard nodes for linear spline
% USAGE
%   x=LINNODE(breaks,evennum)
%
% Standard nodes are the breakpoints
%
% See also: LINBASE, LINDEF, FUNNODE

% Copyright (c) 1997-2000, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function x=linnode(breaks,evennum)
x=breaks;