%% DISCMOMENTS
%
%  Generates means and standard deviations of discrete distribution.
%
%  Usage
%    [xavg,xstd] = discmoments(w,x)
%  Input
%    w    : n.1 discrete probabilities
%    x    : n.k values for k distinct discrete variates
%  Output
%    xavg : 1.k means
%    xstd : 1.k standard deviations

%  Copyright(c) 2014
%   Mario J. Miranda - miranda.4@osu.edu

function [xavg,xstd] = discmoments(w,x)
xavg = w'*x;
xstd = sqrt(w'*x.^2-xavg.^2);